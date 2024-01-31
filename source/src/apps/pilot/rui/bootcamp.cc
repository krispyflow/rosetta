// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <iostream>
#include <basic/options/option.hh>
#include <devel/init.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>

#include <protocols/moves/PyMOLMover.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>

int
main( int argc, char ** argv ) {
	devel::init( argc, argv );
	utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
	if ( filenames.size() > 0 ) {
        
        core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1]);
        
        core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
        core::Real score = sfxn->score( *mypose );
        
        std::cout << "score: " << score << std::endl;
        
		std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
        
        
        
        core::Real uniform_random_number = numeric::random::uniform();
        
        core::Size N = (mypose->size()) / (mypose->total_residue());
        //… code here to pick the index of a random residue in the Pose
        core::Size randres = static_cast< core::Size > ( uniform_random_number * N + 1 );
        

        protocols::moves::MonteCarlo mc1(*mypose, *sfxn, 1);
        
        protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( *mypose, true, 0 );
        the_observer->pymol().apply( *mypose);
         
        int acceptCount = 0;
        int totalItr = 500;
        
        for (int itr = 0; itr < totalItr; ++itr) {
            core::Real pert1 = numeric::random::gaussian();
            //… code here to get a random number
            core::Real pert2 = numeric::random::gaussian();
            //… code here to get another random number
            core::Real orig_phi = mypose->phi( randres );
            core::Real orig_psi = mypose->psi( randres );
            mypose->set_phi( randres, orig_phi + pert1 );
            mypose->set_psi( randres, orig_psi + pert2 );
            
            core::pack::task::PackerTaskOP repack_task =
            core::pack::task::TaskFactory::create_packer_task( *mypose );
            repack_task->restrict_to_repacking();
            core::pack::pack_rotamers( *mypose, *sfxn, repack_task );
            
            core::kinematics::MoveMap mm;
            mm.set_bb( true );
            mm.set_chi( true );

            core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
            core::optimization::AtomTreeMinimizer atm;
            //run minimizer:
            atm.run( *mypose, mm, *sfxn, min_opts );
            
            core::pose::Pose copy_pose = *mypose;
            atm.run( copy_pose, mm, *sfxn, min_opts );
            *mypose = copy_pose;

            mc1.boltzmann(*mypose);
            int modul = itr % 100;
            if (modul == 0) {
                if (mc1.accept() == true) {
                    ++acceptCount;
                }
            }
        }
        core::Real acceptRate = acceptCount / (totalItr/100);
        
	} else {
		std::cout << "You did not provide a PDB file with the -in::file::s option" << std::endl;
		return 1;
	}
}

// bin/bootcamp.default.macosclangdebug -in:file:s /Users/apple/Downloads/rosetta/testFile/1ubq.pdb -database /Users/apple/Downloads/rosetta/rosetta/database -ignore_unrecognized_res
