// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/ContactMapEvaluatorCreator.hh
/// @brief  Header for ContactMapEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/ContactMapEvaluatorCreator.hh>

// Package Headers

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

#include <protocols/simple_filters/ContactMapEvaluator.hh>

#include <core/pose/Pose.hh>

// ObjexxFCL Headers

// Utility headers

#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

//// C++ headers

// due to template function


// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>

#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh>
#endif


static basic::Tracer tr( "protocols.evalution.ContactMapEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

ContactMapEvaluatorCreator::~ContactMapEvaluatorCreator() = default;

void ContactMapEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::contact_map );
	OPT( in::file::native );
}

void ContactMapEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;


	if ( option[ OptionKeys::evaluation::contact_map ] ) {
		core::pose::PoseOP native_pose = nullptr;
		if ( option[ in::file::native ].user() ) {
			native_pose = utility::pointer::make_shared< core::pose::Pose >();
			core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
		}

		if ( !native_pose ) {
			tr.Error << "-evaluation::contact_map must be specified with a native!\n";
		} else {
			core::Real max_dist(12);
			core::Size min_seqsep(5);
			eval.add_evaluation(
				utility::pointer::make_shared< ContactMapEvaluator >( *native_pose, max_dist, min_seqsep )
			);
		}
	}


}

std::string ContactMapEvaluatorCreator::type_name() const {
	return "ContactMapEvaluatorCreator";
}

} //namespace
} //namespace
