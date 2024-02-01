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
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/moves/Mover.hh>

int
main( int argc, char ** argv ) {
    devel::init(argc, argv);
    protocols::moves::MoverOP mover(new protocols::bootcamp::BootCampMover);
    protocols::jd2::JobDistributor::get_instance()->go(mover);
        
}

// bin/bootcamp.default.macosclangdebug -in:file:s /Users/apple/Downloads/rosetta/testFile/1ubq.pdb -database /Users/apple/Downloads/rosetta/rosetta/database -ignore_unrecognized_res
