// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/public/rosettaDNA/rosettaDNA.cc
/// @brief Runs DockDesignParser via the JobDistributor v. 2, provides access to protein-DNA interface modeling code in protocols/dna, and uses a custom PDBJobOutputter to add protein-DNA specific information to output pdb files.

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/viewer/viewers.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>

#include <core/scoring/constraints/Constraint.fwd.hh>
#include <utility/excn/Exceptions.hh>
using namespace core;
using namespace scoring;
using namespace constraints;

using namespace protocols;
using namespace dna;
using namespace moves;
using namespace jd2;
using namespace viewer;

static basic::Tracer TR( "apps.public.rosettaDNA" );

/// @brief this "DummyMover" is employed simply to prevent NULL pointer exceptions in case internal classes try to call Mover methods without checking for pointer validity. The JobDistributor/DockDesignParser will/should reassign this pointer to a real Mover.
class DummyMover : public Mover {
public:
	DummyMover() : Mover( "DummyMover" ) {}
	~DummyMover() override= default;
	void apply( Pose & ) override {
		TR.Fatal << "DummyMover::apply() should never have been called!"
			<< " (JobDistributor/Parser should have replaced DummyMover.)" << std::endl;
		utility_exit_with_message("Unimplemented method called.");
	}
	std::string get_name() const override { return "DummyMover"; }
};

void *
my_main( void * )
{
	MoverOP dummy_mover( new DummyMover );

	// loads dna-specific PDBOutput class into the JobDistributor
	JobDistributor::get_instance()->go( dummy_mover, utility::pointer::make_shared< PDBOutput >() );

	TR << "*********************successful completion**************************" << std::endl;
	return nullptr;
}

int
main( int argc, char * argv [] )
{
	try {
		devel::init( argc, argv );
		viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
}
