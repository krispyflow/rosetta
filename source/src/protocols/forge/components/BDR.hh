// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/components/BDR.hh
/// @brief
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_components_BDR_hh
#define INCLUDED_protocols_forge_components_BDR_hh

// unit headers
#include <protocols/forge/components/BDR.fwd.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.fwd.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/components/VarLengthBuild.fwd.hh>

// project headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>

// utility headers

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace components {


class BDR : public protocols::moves::Mover {


private: // typedefs


	typedef protocols::moves::Mover Super;


public: // typedefs


	typedef std::string String;

	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pack::task::TaskFactoryOP TaskFactoryOP;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;

	typedef protocols::forge::build::BuildInstructionOP BuildInstructionOP;
	typedef protocols::forge::build::BuildManager BuildManager;
	typedef protocols::forge::build::BuildManager::Original2Modified Original2Modified;
	typedef protocols::forge::build::Interval Interval;
	typedef protocols::moves::MoverOP MoverOP;

	/// @brief A pair storing a instruction's original interval and a string
	///  denoting the sequence during design.
	/// @remarks This is only for instructions such as SegmentRebuild or
	///  SegmentInsert.  Non-applicable instructions will not have any data
	///  stored here.
	typedef std::pair< Interval, String > OriginalInterval2DesignString;
	typedef utility::vector1< OriginalInterval2DesignString > DesignInfo;


public: // construct/destruct


	/// @brief default constructor
	BDR();


	/// @brief copy constructor
	BDR( BDR const & rval );


	/// @brief default destructor
	~BDR() override;


private: // disallow assignment


	/// @brief copy assignment
	/// @remarks Mover base class prevents this from working properly...
	BDR & operator =( BDR const & rval );


public: // virtual constructors


	/// @brief clone this object
	MoverOP clone() const override;


	/// @brief create this type of object
	MoverOP fresh_instance() const override;


public: // accessors


	/// @brief use full-mer fragments when building the loop? (default false)
	inline
	bool use_fullmer() const {
		return use_fullmer_;
	}


	/// @brief use sequence biased fragments when building the loop? (default false)
	inline
	bool use_sequence_bias() const {
		return use_sequence_bias_;
	}


	/// @brief the maximum allowed linear chainbreak (default 0.07)
	inline
	Real max_linear_chainbreak() const {
		return max_linear_chainbreak_;
	}


	/// @brief the loop mover string to use during centroid build
	///  (default "RemodelLoopMover")
	/// @remarks set to either a string the create_loop_mover() LoopMoverFactory
	///  recognizes or the "RemodelLoopMover"
	inline
	String const & centroid_loop_mover_str() const {
		return centroid_loop_mover_str_;
	}


	/// @brief redesign the neighborhood around the loop?  if false, then just
	///  repacks during the design phase (default true)
	inline
	bool redesign_loop_neighborhood() const {
		return redesign_loop_neighborhood_;
	}


	/// @brief name of the resfile to use during design-refine; empty string
	///  implies no resfile
	inline
	String const & resfile() const {
		return resfile_;
	}


	/// @brief the number of design-refine cycles to perform, default 3
	inline
	core::Size dr_cycles() const {
		return dr_cycles_;
	}


	/// @brief the centroid level score function, default "remodel_cen"
	ScoreFunction const & centroid_scorefunction() const;


	/// @brief the full-atom level score function, default score12
	ScoreFunction const & fullatom_scorefunction() const;


public: // mutators


	/// @brief add instruction to the manager of this BDR (no copy)
	/// @param[in] bi BuildInstruction
	/// @param[in] aa_during_design_refine The allowed amino acid sequence
	///  during design.  Only applicable to BuildInstructions like
	///  SegmentRebuild and SegmentInsert.  Make sure the length of this
	///  string matches up properly, and remember to use any special characters,
	///  e.g. the insertion character for SegmentInsert
	void add_instruction(
		BuildInstructionOP bi,
		String const & aa_during_design_refine = String()
	);


	/// @brief create directed dependency between two instructions
	void create_directed_dependency(
		BuildInstructionOP u,
		BuildInstructionOP v
	);


	/// @brief use full-mer fragments when building the loop?
	inline
	void use_fullmer( bool const flag ) {
		use_fullmer_ = flag;
	}


	/// @brief use sequence biased fragments when building the loop? (default false)
	inline
	void use_sequence_bias( bool const flag ) {
		use_sequence_bias_ = flag;
	}


	/// @brief the maximum allowed linear chainbreak
	inline
	void max_linear_chainbreak( Real const threshold ) {
		max_linear_chainbreak_ = threshold;
	}


	/// @brief the loop mover string to use during centroid build
	/// @remarks set to either a string the create_loop_mover() LoopMoverFactory
	///  recognizes or the "RemodelLoopMover"
	inline
	void centroid_loop_mover_str( String const & loop_mover_str ) {
		centroid_loop_mover_str_ = loop_mover_str;
	}


	/// @brief redesign the neighborhood around the loop?  if false, then just
	///  repacks during the design phase
	inline
	void redesign_loop_neighborhood( bool const flag ) {
		redesign_loop_neighborhood_ = flag;
	}


	/// @brief name of the resfile to use during design-refine; empty string
	///  implies no resfile
	inline
	void resfile( String const & filename ) {
		resfile_ = filename;
	}

	/// @brief set the number of design-refine cycles to perform
	/// @remarks set this to 0 to skip design-refine
	inline
	void dr_cycles( core::Size const cycles ) {
		dr_cycles_ = cycles;
	}


	/// @brief set the centroid level score function
	void centroid_scorefunction( ScoreFunction const & sfx );


	/// @brief set the centroid level score function
	void centroid_scorefunction( ScoreFunctionOP sfx );


	/// @brief set the full-atom level score function
	void fullatom_scorefunction( ScoreFunction const & sfx );


	/// @brief set the full-atom level score function
	void fullatom_scorefunction( ScoreFunctionOP sfx );


public: // calculators


	/// @brief the name for the loops' neighborhood calculator
	inline
	static
	String neighborhood_calc_name() {
		return "BDR_loops_neighborhood_calc";
	}


	/// @brief the name for the loops' buried unsatisfied polars
	///  calculator
	inline
	static
	String loops_buns_polar_calc_name() {
		return "BDR_loops_bup";
	}


	/// @brief the name for the loops' neighborhood buried unsatisfied polars
	///  calculator
	inline
	static
	String neighborhood_buns_polar_calc_name() {
		return "BDR_loops_neighborhood_bup";
	}


public: // virtual main methods


	/// @brief apply defined moves to given Pose
	void apply( Pose & pose ) override;

	std::string get_name() const override;

private: // protocol methods


	/// @brief run the centroid level build stage
	/// @return true if loop closed, false otherwise
	bool centroid_build(
		Pose & pose
	);


	/// @brief run the design-refine stage
	/// @return currently always true
	bool design_refine(
		Pose & pose
	);


	/// @brief return a TaskFactory useable as a starting point for either
	///  design or refinement
	TaskFactoryOP generic_taskfactory();


private: // design methods


	/// @brief process a continuous design string, adding appropriate operations
	///  to the TaskFactory
	void process_continuous_design_string(
		Interval const & original_interval,
		String const & design_str,
		Original2Modified const & original2modified_interval_endpoints,
		TaskFactoryOP design_tf
	);


	/// @brief process a design string containing an insert, adding appropriate
	///  operations to the TaskFactory
	void process_insert_design_string(
		Interval const & original_interval,
		String const & design_str,
		Original2Modified const & original2modified_interval_endpoints,
		TaskFactoryOP design_tf
	);


	/// @brief return a boolean vector specifying allowed a.a. when designing
	///  on the surface
	static
	utility::vector1< bool > const & allowed_surface_aa();


private: // data


	/// @brief the BuildManager
	BuildManager manager_;


	/// @brief Stores instructions' original interval and a string denoting
	///  the sequence during design.
	/// @remarks This is only for instructions such as SegmentRebuild or
	///  SegmentInsert.  Non-applicable instructions will not have any data
	///  stored here.
	DesignInfo design_info_;


	/// @brief use full-mer fragments when building the loop? (default false)
	bool use_fullmer_;


	/// @brief use sequence biased fragments when building the loop? (default false)
	bool use_sequence_bias_;


	/// @brief the maximum allowed linear chainbreak (default 0.07)
	Real max_linear_chainbreak_;


	/// @brief the loop mover string to use during centroid build
	///  (default "RemodelLoopMover")
	/// @remarks set to either a string the create_loop_mover() LoopMoverFactory
	///  recognizes or the "RemodelLoopMover"
	String centroid_loop_mover_str_;


	/// @brief redesign the neighborhood around the loop?  if false, then just
	///  repacks during the design phase (default true)
	bool redesign_loop_neighborhood_;


	/// @brief name of the resfile to use during design-refine; empty string
	///  implies no resfile
	String resfile_;


	/// @brief the number of design-refine cycles to perform, default 3
	core::Size dr_cycles_;


	/// @brief the centroid scorefunction to use, default "remodel_cen"
	ScoreFunctionOP centroid_sfx_;


	/// @brief the full-atom scorefunction to use, default score12
	ScoreFunctionOP fullatom_sfx_;


private: // per-stage movers


	/// @brief VLB for centroid_build
	/// @remarks Store it here instead of re-instantiation each centroid_build()
	///  call so fragment caching works.
	/// @warning Remember to set this to null if the BuildManager changes.
	VarLengthBuildOP vlb_;


};


} // namespace components
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_components_BDR_HH */
