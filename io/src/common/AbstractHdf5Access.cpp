/*

Copyright (c) 2005-2023, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "AbstractHdf5Access.hpp"
#include "Exception.hpp"

bool AbstractHdf5Access::DoesDatasetExist(const std::string& rDatasetName)
{
    // This is a nice method for testing existence, introduced in HDF5 1.8.0
    htri_t dataset_status = H5Lexists(mFileId, rDatasetName.c_str(), H5P_DEFAULT);
    return (dataset_status > 0);
}

void AbstractHdf5Access::SetUnlimitedDatasetId()
{
    // Now deal with the unlimited dimension

    // In terms of an Unlimited dimension dataset:
    // * Files pre - r16738 (inc. Release 3.1 and earlier) use simply "Time" for "Data"'s unlimited variable.
    // * Files generated by r16738 - r18257 used "<DatasetName>_Time" for "<DatasetName>"'s unlimited variable,
    //   - These are not to be used and there is no backwards compatibility for them, since they weren't in a release.
    // * Files post r18257 (inc. Release 3.2 onwards) use "<DatasetName>_Unlimited" for "<DatasetName>"'s
    //   unlimited variable,
    //   - a new attribute "Name" has been added to the Unlimited Dataset to allow it to assign
    //     any name to the unlimited variable. Which can then be easily read by Hdf5DataReader.
    //   - if this dataset is missing we look for simply "Time" to remain backwards compatible with Releases <= 3.1

    if (DoesDatasetExist(mDatasetName + "_Unlimited"))
    {
        mUnlimitedDatasetId = H5Dopen(mFileId, (mDatasetName + "_Unlimited").c_str(), H5P_DEFAULT);
        hid_t name_attribute_id = H5Aopen_name(mUnlimitedDatasetId, "Name");
        hid_t unit_attribute_id = H5Aopen_name(mUnlimitedDatasetId, "Unit");

        hid_t attribute_type = H5Aget_type(name_attribute_id);

        // Read into it.
        char* string_array = (char*)malloc(sizeof(char) * MAX_STRING_SIZE);
        H5Aread(name_attribute_id, attribute_type, string_array);
        std::string name_string(&string_array[0]);
        mUnlimitedDimensionName = name_string;

        H5Aread(unit_attribute_id, attribute_type, string_array);
        std::string unit_string(&string_array[0]);
        mUnlimitedDimensionUnit = unit_string;

        free(string_array);
        H5Tclose(attribute_type);
        H5Aclose(name_attribute_id);
        H5Aclose(unit_attribute_id);
    }
    else if (DoesDatasetExist("Time"))
    {
        mUnlimitedDimensionName = "Time";
        mUnlimitedDimensionUnit = "msec";
        mUnlimitedDatasetId = H5Dopen(mFileId, mUnlimitedDimensionName.c_str(), H5P_DEFAULT);
    }
    else
    {
        NEVER_REACHED;
    }
    mIsUnlimitedDimensionSet = true;
}

AbstractHdf5Access::AbstractHdf5Access(const std::string& rDirectory,
                                       const std::string& rBaseName,
                                       const std::string& rDatasetName,
                                       bool makeAbsolute)
        : mBaseName(rBaseName),
          mDatasetName(rDatasetName),
          mIsDataComplete(true),
          mIsUnlimitedDimensionSet(false)
{
    RelativeTo::Value relative_to;
    if (makeAbsolute)
    {
        relative_to = RelativeTo::ChasteTestOutput;
    }
    else
    {
        relative_to = RelativeTo::Absolute;
    }
    mDirectory.SetPath(rDirectory, relative_to);
}

AbstractHdf5Access::AbstractHdf5Access(const FileFinder& rDirectory,
                                       const std::string& rBaseName,
                                       const std::string& rDatasetName)
        : mBaseName(rBaseName),
          mDatasetName(rDatasetName),
          mDirectory(rDirectory),
          mIsDataComplete(true),
          mIsUnlimitedDimensionSet(false)
{
}

AbstractHdf5Access::~AbstractHdf5Access()
{
}

bool AbstractHdf5Access::IsDataComplete()
{
    return mIsDataComplete;
}

std::vector<unsigned> AbstractHdf5Access::GetIncompleteNodeMap()
{
    return mIncompleteNodeIndices;
}

std::string AbstractHdf5Access::GetUnlimitedDimensionName()
{
    return mUnlimitedDimensionName;
}

std::string AbstractHdf5Access::GetUnlimitedDimensionUnit()
{
    return mUnlimitedDimensionUnit;
}

void AbstractHdf5Access::SetMainDatasetRawChunkCache()
{
    // 128 M cache for raw data. 12799 is a prime number which is 100 x larger
    // than the number of 1 MB chunks (the largest we have tested with) which
    // could fit into the 128 M cache. 0.75 for the last argument is the default.
    // See: http://www.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-SetChunkCache

    hsize_t max_objects_in_chunk_cache = 12799u;
    hsize_t max_bytes_in_cache = 128u * 1024u * 1024u;
#if H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 8 && H5_VERS_RELEASE >= 3 // HDF5 1.8.3+
    // These methods set the cache on a dataset basis
    hid_t dapl_id = H5Dget_access_plist(mVariablesDatasetId);
    H5Pset_chunk_cache(dapl_id,
                       max_objects_in_chunk_cache,
                       max_bytes_in_cache,
                       H5D_CHUNK_CACHE_W0_DEFAULT);
#else
    // These older methods set the cache on a file basis
    hid_t fapl_id = H5Fget_access_plist(mFileId);
    H5Pset_cache(fapl_id,
                 0, // unused
                 max_objects_in_chunk_cache,
                 max_bytes_in_cache,
                 0.75); // default value
#endif
}
