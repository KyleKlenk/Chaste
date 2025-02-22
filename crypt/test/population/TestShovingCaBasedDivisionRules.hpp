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

#ifndef TESTSHOVINGCABASEDDIVISIONRULES_HPP_
#define TESTSHOVINGCABASEDDIVISIONRULES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "SmartPointers.hpp"
#include "ArchiveOpener.hpp"
#include "CellsGenerator.hpp"
#include "CaBasedCellPopulation.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCaBasedDivisionRule.hpp"
#include "ShovingCaBasedDivisionRule.hpp"
#include "CryptShovingCaBasedDivisionRule.hpp"
#include "PottsMeshGenerator.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * \todo Merge this test suite into cell_based/test/population/TestCaBasedDivisionRules once
 * the contents of notforrelease_cell_based/src/population/division_rules are moved into
 * cell_based/src/population/division_rules or crypt/src/population (#2731)
 */
class TestShovingCaBasedDivisionRules : public AbstractCellBasedTestSuite
{
public:

    /**
     * In this test we create a new ShovingCaBasedDivisionRule, divide a cell with it
     * and check that the new cells are in the correct locations. First, we test where
     * there is space around the cells. This is the default setup.
     */
    void TestAddCellWithShovingBasedDivisionRule()
    {
        // Create a simple Potts mesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create 9 cells in the central nodes
        std::vector<unsigned> location_indices;
        for (unsigned row=1; row<4; row++)
        {
            location_indices.push_back(1+row*5);
            location_indices.push_back(2+row*5);
            location_indices.push_back(3+row*5);
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Check the cell locations
        unsigned cell_locations[9] = {6, 7, 8, 11, 12, 13, 16, 17, 18};
        unsigned index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter),cell_locations[index])
            ++index;
        }

        // Make a new cell to add
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetCellProliferativeType(p_stem_type);
        p_new_cell->SetBirthTime(-1);

        // Set the division rule for our population to be the shoving division rule
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule_to_set(new ShovingCaBasedDivisionRule<2>());
        cell_population.SetCaBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population and try to add new cell by dividing cell at site 0
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule = cell_population.GetCaBasedDivisionRule();

        // Select central cell
        CellPtr p_cell_12 = cell_population.GetCellUsingLocationIndex(12);

        // The ShovingCaBasedDivisionRule method IsRoomToDivide() always returns true
        TS_ASSERT_EQUALS((p_division_rule->IsRoomToDivide(p_cell_12, cell_population)), true);

        /*
         * Test adding the new cell to the population; this calls CalculateDaughterNodeIndex().
         * The new cell moves into node 13.
         */
        cell_population.AddCell(p_new_cell, p_cell_12);

        // Now check the cells are in the correct place
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 10u);

        // Note the cell originally on node 13 has been shoved to node 14 and the new cell is on node 13
        unsigned new_cell_locations[10] = {6, 7, 8, 11, 12, 14, 16, 17, 18, 13};
        index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), new_cell_locations[index])
            ++index;
        }
    }

    /**
     * In this test of ShovingCaBasedDivisionRule we check the case where there is
     * no room to divide without the cells being shoved to the edge of the mesh.
     */
    void TestAddCellWithShovingBasedDivisionRuleAndShovingRequired()
    {
        // Create a simple Potts mesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create 25 cells, one for each node
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<25; index++)
        {
            location_indices.push_back(index);
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Make a new cell to add
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetCellProliferativeType(p_stem_type);
        p_new_cell->SetBirthTime(-1);

        // Set the division rule for our population to be the shoving division rule
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule_to_set(new ShovingCaBasedDivisionRule<2>());
        cell_population.SetCaBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population and try to add new cell by dividing cell at site 0;
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule = cell_population.GetCaBasedDivisionRule();

        // Select central cell
        CellPtr p_cell_12 = cell_population.GetCellUsingLocationIndex(12);

        // Try to divide but cant as hit boundary
        TS_ASSERT_THROWS_THIS(p_division_rule->CalculateDaughterNodeIndex(p_new_cell, p_cell_12, cell_population),
            "Cells reaching the boundary of the domain. Make the Potts mesh larger.");
    }

    void TestArchivingShovingCaBasedDivisionRule()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ShovingCaBasedDivisionRule.arch";

        {
            ShovingCaBasedDivisionRule<2> division_rule;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer to most abstract class possible
            AbstractCaBasedDivisionRule<2>* const p_division_rule = &division_rule;
            output_arch << p_division_rule;
            ofs.close();
        }

        {
            AbstractCaBasedDivisionRule<2>* p_division_rule;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_division_rule;
            ifs.close();

            TS_ASSERT(p_division_rule != NULL);

            // Tidy up
            delete p_division_rule;
        }
    }

    /**
     * In this test we create a new CryptShovingCaBasedDivisionRule, divide a cell with it
     * and check that the new cells are in the correct locations. First, we test where
     * there is space around the cells. This is the default setup.
     */
    void TestAddCellWithCryptShovingBasedDivisionRule()
    {
        // Create a simple Potts mesh
        PottsMeshGenerator<2> generator(3, 0, 0, 4, 0, 0,1,0,0,false, true); // Periodic in x
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create 6 cells in the bottom 2 rows
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<6; index++)
        {
            location_indices.push_back(index);
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size()); // Note all cells are stem cells by default.

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Check the cell locations
        unsigned cell_locations[6] = {0,1,2,3,4,5};
        unsigned index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter),cell_locations[index])
            ++index;
        }

        // Make a new cell to add
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetCellProliferativeType(p_transit_type);
        p_new_cell->SetBirthTime(-1);

        // Set the division rule for our population to be the cryot shoving division rule
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule_to_set(new CryptShovingCaBasedDivisionRule());
        cell_population.SetCaBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population and try to add new cell by dividing cell at site 0;
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule = cell_population.GetCaBasedDivisionRule();

        // Select left cell in bottom row
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);

        // The CryptShovingCaBasedDivisionRule method IsRoomToDivide() always returns true
        TS_ASSERT((p_division_rule->IsRoomToDivide(p_cell_0,cell_population)));

        /*
         * Test adding the new cell in the population; this calls CalculateDaughterNodeIndex().
         * The new cell moves into node 3. This is because stem cells always divide upwards
         */
        cell_population.AddCell(p_new_cell, p_cell_0);

        // Now check the cells are in the correct place
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 7u);

        // Note the cell on node 3 has been shoved to node 6 and the new cell is on node 3
        unsigned new_cell_locations[7] = {0,1,2,6,4,5,3};
        index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter),new_cell_locations[index])
            ++index;
        }

        //Now divide a non stem cell
        // Select transit cell we just added
        CellPtr p_cell_3 = cell_population.GetCellUsingLocationIndex(3);

        FixedG1GenerationalCellCycleModel* p_model_2 = new FixedG1GenerationalCellCycleModel();
        CellPtr p_new_cell_2(new Cell(p_state, p_model_2));
        p_new_cell_2->SetCellProliferativeType(p_transit_type);
        p_new_cell_2->SetBirthTime(-1);

        /*
         * Test adding the new cell in the population; this calls CalculateDaughterNodeIndex().
         * The new cell moves into node 7.
         */
        cell_population.AddCell(p_new_cell_2, p_cell_3);

        // Now check the cells are in the correct place
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);

        // Note the cell on node 4 has been shoved to node 7 and the new cell is on node 4
        unsigned new_cell_locations_2[8] = {0,1,2,6,7,5,3,4};
        index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter),new_cell_locations_2[index])
            ++index;
        }
    }

    /**
     * In this test of CryptShovingCaBasedDivisionRule we check the case where there is
     * no room to divide without the cells being shoved to the edge of the mesh.
     */
    void TestAddCellWithCryptShovingBasedDivisionRuleAndShovingRequired()
    {
        // Create a simple Potts mesh
        PottsMeshGenerator<2> generator(3, 0, 0, 3, 0, 0, 1, 0, 0, false, true); // x periodic
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create 9 cells, one for each node
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<9; index++)
        {
            location_indices.push_back(index);
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Make a new cell to add
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetCellProliferativeType(p_stem_type);
        p_new_cell->SetBirthTime(-1);

        // Set the division rule for our population to be the crypt shoving division rule
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule_to_set(new CryptShovingCaBasedDivisionRule());
        cell_population.SetCaBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population and try to add new cell by dividing cell at site 0
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule = cell_population.GetCaBasedDivisionRule();

        // Select bottom left cell
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);

        // Can't divide without shoving into top
        TS_ASSERT_THROWS_THIS(p_division_rule->CalculateDaughterNodeIndex(p_new_cell, p_cell_0, cell_population),
            "Cells reaching the top of the crypt need to increase length to at least double the sloughing height.");
    }

    void TestArchivingCryptShovingCaBasedDivisionRule()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CryptShovingCaBasedDivisionRule.arch";

        {
            CryptShovingCaBasedDivisionRule division_rule;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer to most abstract class possible
            AbstractCaBasedDivisionRule<2>* const p_division_rule = &division_rule;
            output_arch << p_division_rule;
            ofs.close();
        }

        {
            AbstractCaBasedDivisionRule<2>* p_division_rule;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_division_rule;
            ifs.close();

            TS_ASSERT(p_division_rule != NULL);

            // Tidy up
            delete p_division_rule;
        }
    }
};

#endif /*TESTSHOVINGCABASEDDIVISIONRULES_HPP_*/
