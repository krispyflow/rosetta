// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author shilei

#include <iostream>

#include <protocols/constraint_movers/ConstraintSetMover.hh>

#include <devel/init.hh>

//pose
#include <core/pose/Pose.hh>

//score
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//job distribution

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

//Options
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>


//Auto Headers
#include <core/import_pose/import_pose.hh>

//docking
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/metrics.hh>

//jumps
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Conformation.hh>

//matrix
#include <numeric/xyzMatrix.hh>

//id


//silent

//options
#include <basic/options/option_macros.hh>


using namespace core;
using namespace core::scoring;
using namespace std;
using utility::vector1;
using basic::options::option;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace pose;
using std::string;

///work/shilei/rosetta/rosetta_source/src/apps/pilot/wendao/test_bbmc.cc
OPT_1GRP_KEY(String,read_pose_jump_orientation_repack,outtag)
OPT_1GRP_KEY(RealVector,read_pose_jump_orientation_repack,jump_orientation)


void transform_pose( core::pose::Pose & pose,utility::vector1<core::Real> const & t )
{
	using core::Size;
	using core::Real;
	using namespace core::id;
	using namespace core::scoring;
	using core::kinematics::Jump;

	//add rb changes
	numeric::xyzMatrix<core::Real> rotation_matrix=numeric::xyzMatrix<Real>::rows(t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9]);
	numeric::xyzVector<core::Real> translation_vector=numeric::xyzVector<core::Real>(t[10],t[11],t[12]);
	core::kinematics::Jump tmpJump=pose.jump(1);
	tmpJump.set_rotation(rotation_matrix);
	tmpJump.set_translation(translation_vector);
	pose.conformation().set_jump(1,tmpJump);
	//pose.conformation().set_jump_now(1,tmpJump);

}

///////////////////////////////////////////////////////////////////////////////
void run() {

	//read input pose on the every node, rather than communicating
	Pose pose,ref_pose;
	std::string pdbname;
	if ( basic::options::option[ basic::options::OptionKeys::in::file::s ].user() ) {
		pdbname=basic::options::option[ basic::options::OptionKeys::in::file::s ]()[1];
		core::import_pose::pose_from_file( pose, pdbname.c_str() , core::import_pose::PDB_file);
		ref_pose=pose;
	} else {
		throw ( CREATE_EXCEPTION(utility::excn::BadInput, "expected -s for this app") );
	}

	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		core::import_pose::pose_from_file( ref_pose, basic::options::option[ basic::options::OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
	}

	utility::vector1< core::Real > saved_transformations_;
	if ( !basic::options::option[ basic::options::OptionKeys::read_pose_jump_orientation_repack::jump_orientation].user() ||
			basic::options::option[ basic::options::OptionKeys::read_pose_jump_orientation_repack::jump_orientation].size()!=12 ) {
		throw ( CREATE_EXCEPTION(utility::excn::BadInput, "expected jump_orientation be provided as a vector of 12") );
	} else {
		for ( core::Size i=1; i <= 12; i++ ) {
			saved_transformations_.push_back(basic::options::option[ basic::options::OptionKeys::read_pose_jump_orientation_repack::jump_orientation][i]);
		}
	}

	//for (core::Size i=1; i <= 12; i++) {
	// std::cout << saved_transformations_[i] << std::endl;
	//}
	//
	//set up fold-tree
	std::cout << "foldtree 0: "<< pose.fold_tree() << endl;
	utility::vector1< core::Size > movable_jumps_ = utility::tools::make_vector1<core::Size>(1);
	//protocols::docking::DockJumps movable_jumps_(1);
	//protocols::docking::setup_foldtree( pose, utility::to_string("_"), movable_jumps_);
	//std::cout << "foldtree 1: "<< pose.fold_tree() << endl;

	//pose.dump_pdb( "original.pdb");
	//std::cout << " original: ";
	//std::cout << pose.jump(1).get_rotation()(1,1) << " ";
	//std::cout << pose.jump(1).get_rotation()(1,2) << " ";
	//std::cout << pose.jump(1).get_rotation()(1,3) << " ";
	//std::cout << pose.jump(1).get_rotation()(2,1) << " ";
	//std::cout << pose.jump(1).get_rotation()(2,2) << " ";
	//std::cout << pose.jump(1).get_rotation()(2,3) << " ";
	//std::cout << pose.jump(1).get_rotation()(3,1) << " ";
	//std::cout << pose.jump(1).get_rotation()(3,2) << " ";
	//std::cout << pose.jump(1).get_rotation()(3,3) << " ";
	//std::cout << pose.jump(1).get_translation()(1) << " ";
	//std::cout << pose.jump(1).get_translation()(2) << " ";
	//std::cout << pose.jump(1).get_translation()(3) << std::endl;
	transform_pose(pose,saved_transformations_);
	Real rms = core::scoring::CA_rmsd(pose,ref_pose);
	//pose.dump_pdb( "transform.pdb");
	std::cout << "rms_after_transform: " << rms << " ";
	//std::cout << " transform: ";
	//std::cout << pose.jump(1).get_rotation()(1,1) << " ";
	//std::cout << pose.jump(1).get_rotation()(1,2) << " ";
	//std::cout << pose.jump(1).get_rotation()(1,3) << " ";
	//std::cout << pose.jump(1).get_rotation()(2,1) << " ";
	//std::cout << pose.jump(1).get_rotation()(2,2) << " ";
	//std::cout << pose.jump(1).get_rotation()(2,3) << " ";
	//std::cout << pose.jump(1).get_rotation()(3,1) << " ";
	//std::cout << pose.jump(1).get_rotation()(3,2) << " ";
	//std::cout << pose.jump(1).get_rotation()(3,3) << " ";
	//std::cout << pose.jump(1).get_translation()(1) << " ";
	//std::cout << pose.jump(1).get_translation()(2) << " ";
	//std::cout << pose.jump(1).get_translation()(3) << std::endl;

	//repack interface and computer Isc, rms
	core::scoring::ScoreFunctionOP scorefxn_cen( ScoreFunctionFactory::create_score_function("interchain_cen") );
	core::scoring::ScoreFunctionOP docking_scorefxn_high_ = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	core::scoring::ScoreFunctionOP cst_score_( core::scoring::ScoreFunctionFactory::create_score_function("empty") );

	//add cst to the scoring
	if ( basic::options::option[basic::options::OptionKeys::constraints::cst_file].user() ) {
		protocols::constraint_movers::ConstraintSetMoverOP docking_constraint_( new protocols::constraint_movers::ConstraintSetMover() );
		Real cst_weight_=basic::options::option[basic::options::OptionKeys::constraints::cst_weight];
		docking_constraint_->apply(pose);
		cst_score_->set_weight( core::scoring::atom_pair_constraint, cst_weight_ );
	}

	// local_tf = new TaskFactory( );
	// local_tf->push_back( new RestrictToRepacking );
	// local_tf->push_back( new InitializeFromCommandline );
	// local_tf->push_back( new IncludeCurrent );
	// local_tf->push_back( new NoRepackDisulfides );
	//      core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot = new core::pack::rotamer_set::UnboundRotamersOperation();
	//        unboundrot->initialize_from_command_line();
	//        operation::AppendRotamerSetOP unboundrot_operation = new operation::AppendRotamerSet( unboundrot );
	//        local_tf->push_back( unboundrot_operation );
	//        core::pack::dunbrack::load_unboundrot(pose);

	//protocols::docking::DockingHighResLegacyOP docking_highres_mover_ = new protocols::docking::DockingHighResLegacy( movable_jumps_, docking_scorefxn_high_, docking_scorefxn_high_);
	protocols::docking::DockingProtocolOP docking_highres_mover_( new protocols::docking::DockingProtocol( movable_jumps_, false, true, true, scorefxn_cen, docking_scorefxn_high_) );

	//docking_highres_mover_->set_trans_magnitude(0.00001);
	//docking_highres_mover_->set_rot_magnitude(0.00001);
	//std::cout << "trans: "<< docking_highres_mover_->get_trans_magnitude() << std::endl;
	//std::cout << "rot: "<< docking_highres_mover_->get_rot_magnitude() << std::endl;
	// std::cout << " before docking: ";
	// std::cout << pose.jump(1).get_rotation()(1,1) << " ";
	// std::cout << pose.jump(1).get_rotation()(1,2) << " ";
	// std::cout << pose.jump(1).get_rotation()(1,3) << " ";
	// std::cout << pose.jump(1).get_rotation()(2,1) << " ";
	// std::cout << pose.jump(1).get_rotation()(2,2) << " ";
	// std::cout << pose.jump(1).get_rotation()(2,3) << " ";
	// std::cout << pose.jump(1).get_rotation()(3,1) << " ";
	// std::cout << pose.jump(1).get_rotation()(3,2) << " ";
	// std::cout << pose.jump(1).get_rotation()(3,3) << " ";
	// std::cout << pose.jump(1).get_translation()(1) << " ";
	// std::cout << pose.jump(1).get_translation()(2) << " ";
	// std::cout << pose.jump(1).get_translation()(3) << std::endl;
	docking_highres_mover_->apply(pose);
	// std::cout << " after docking: ";
	// std::cout << pose.jump(1).get_rotation()(1,1) << " ";
	// std::cout << pose.jump(1).get_rotation()(1,2) << " ";
	// std::cout << pose.jump(1).get_rotation()(1,3) << " ";
	// std::cout << pose.jump(1).get_rotation()(2,1) << " ";
	// std::cout << pose.jump(1).get_rotation()(2,2) << " ";
	// std::cout << pose.jump(1).get_rotation()(2,3) << " ";
	// std::cout << pose.jump(1).get_rotation()(3,1) << " ";
	// std::cout << pose.jump(1).get_rotation()(3,2) << " ";
	// std::cout << pose.jump(1).get_rotation()(3,3) << " ";
	// std::cout << pose.jump(1).get_translation()(1) << " ";
	// std::cout << pose.jump(1).get_translation()(2) << " ";
	// std::cout << pose.jump(1).get_translation()(3) << std::endl;
	// std::cout << "foldtree 3: "<< pose.fold_tree() << endl;

	//compute the energy
	double CstScore=0.0;
	if ( basic::options::option[basic::options::OptionKeys::constraints::cst_file].user() ) {
		(*cst_score_)(pose);
		CstScore= pose.energies().total_energies().dot( cst_score_->weights() );
	} else {
		CstScore=-999;
	}

	(*docking_scorefxn_high_)(pose);
	double totalScore= pose.energies().total_energies().dot( docking_scorefxn_high_->weights() );
	double Isc=protocols::docking::calc_interaction_energy(pose,docking_scorefxn_high_,movable_jumps_);
	rms = core::scoring::CA_rmsd(pose,ref_pose);
	std::cout << "rms_after_repack: " << rms << " Isc: " << Isc << " totalScore: " << totalScore << " CstScore: " << CstScore << std::endl;

	std::string outtag= basic::options::option[basic::options::OptionKeys::read_pose_jump_orientation_repack::outtag];
	std::string pdbfile=outtag+"transform_repack.pdb";
	pose.dump_pdb( pdbfile );
}

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char * argv [] ) {
	try {
		// initialize option and random number system
		NEW_OPT(read_pose_jump_orientation_repack::outtag, "output pdb tag","");
		NEW_OPT(read_pose_jump_orientation_repack::jump_orientation, "Jump orientation",utility::vector1<core::Real>());
		//NEW_OPT come before devel::init

		devel::init( argc, argv );

		run();

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
}
