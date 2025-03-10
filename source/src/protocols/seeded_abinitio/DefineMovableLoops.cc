// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/seeded_abinitio/
/// @brief specific to seeded abinitio protocol
///  this class finds loops after seeded abinitio that should be closed
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#include <protocols/seeded_abinitio/DefineMovableLoops.hh>
#include <protocols/seeded_abinitio/DefineMovableLoopsCreator.hh>
#include <protocols/seeded_abinitio/SeededAbinitio_util.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/variant_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>

//protocols
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/chemical/VariantType.hh>

// C++ headers
#include <string>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>

//parser
#include <utility/tag/Tag.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.hh>

//util
#include <utility/vector1.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using namespace core;
using namespace protocols::seeded_abinitio;
static basic::Tracer TR( "protocols.seeded_abinitio.DefineMovableLoops" );


namespace protocols {
namespace seeded_abinitio {

using namespace protocols::moves;
using namespace core;




DefineMovableLoops::~DefineMovableLoops() = default;

DefineMovableLoops::DefineMovableLoops() :
	protocols::moves::Mover( DefineMovableLoops::mover_name() )
{
	use_cutpoints_ = true;
	secstructure_ = "";
	chains_.clear();
	chainbreakweights_ = true;
	secstructure_ = "";
	chains_.clear();
	seed_vector_.clear();
}

bool
DefineMovableLoops::chainbreakweights(){
	return chainbreakweights_;
}

void
DefineMovableLoops::add_chainbreakweights( bool acbw ){
	chainbreakweights_ = acbw;
}

void
DefineMovableLoops::use_cutpoints( bool uc ){
	use_cutpoints_ = uc;
}

bool
DefineMovableLoops::use_cutpoints(){
	return use_cutpoints_;
}


protocols::loops::LoopsOP
DefineMovableLoops::find_loops(   pose::Pose & pose,
	std::string secstruct,
	core::Size offset,//first position to start from
	protocols::loops::Loops seeds // change back to OP eventually...
){

	//bool use_seeds( seeds );
	bool use_seeds = seeds.size() > 0;
	core::Size end = offset + secstruct.length() - 1;//double check....
	utility::vector1< core::Size > adjusted_cutpoints;

	//curate the cutpoint list from the cutpoints between the chain ends
	if ( use_cutpoints() ) {
		core::kinematics::FoldTreeOP ft( new kinematics::FoldTree( pose.fold_tree() ) );
		utility::vector1< core::Size > cutpoints =  ft->cutpoints();
		utility::vector1< core::Size > chain_ends = pose.conformation().chain_endings();
		//adding last chains' end too, since it isnt included in chain_endings( and to avoid seg fault beloew)
		chain_ends.push_back( pose.size() );

		//debug
		for ( core::Size i = 1; i <= chain_ends.size(); ++i ) {
			TR.Debug<<"chain endings: " << chain_ends[i] <<std::endl;
		}

		//adding together relevant cutpoints
		for ( core::Size cut_it = 1; cut_it <= cutpoints.size(); ++ cut_it ) {
			TR.Debug <<"cutpoints of current foldtree: "<< cutpoints[cut_it] <<std::endl;

			for ( core::Size ends_it = 1; ends_it <= chain_ends.size(); ++ends_it ) {
				if ( (cutpoints[cut_it] != chain_ends[ends_it])   &&  (cutpoints[cut_it] <= end )  &&  ( cutpoints[cut_it] >= offset) ) {
					adjusted_cutpoints.push_back( cutpoints[cut_it] );
					TR <<"adjusted cutpoint "<< cutpoints[cut_it]  << std::endl;
				}
			}
		}
	}

	//put loop regions together then pull the ones out that contain a cutpoint if use_cutpoints is specified
	protocols::loops::Loops found_loops;

	TR.Debug <<"sec. strc: "<< secstruct <<std::endl;
	//char ss;
	utility::vector1 < core::Size > individual_loop;
	core::Size cut = 0;

	for ( core::Size ss_i = 0; ss_i < secstruct.length() ; ++ss_i ) {
		char ss = secstruct[ss_i];

		if ( ss == 'L' ) {

			//is there a seed and if so, only add to loop if the current indx is not part of it
			if ( use_seeds && !seeds.is_loop_residue( ss_i + offset ) ) {
				individual_loop.push_back( ss_i + offset );
			}

			if ( !use_seeds ) {
				individual_loop.push_back( ss_i + offset );
			}
			TR.Debug <<"use cutpoint " << use_cutpoints_ << " cut: " << cut <<" iterator with offset: " << ss_i + offset <<std::endl;
			if ( use_cutpoints_ && is_cut( adjusted_cutpoints, ss_i + offset ) ) {
				cut = ss_i + offset;
				TR  <<"cut: " << cut << std::endl;
			}
		}//end if statemnt is loop


		if ( ss != 'L' ) {

			//for chainbreak weights ///// this is only added if use cutpoints is set to true!!!!
			if ( chainbreakweights() && cut > 0 ) {
				TR <<"adding chainbreak type"<< std::endl;
				core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cut );
				core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cut+1 );
			}
			//take loop container and add as loop, reset temporariy vector with loop positions
			if ( individual_loop.size() > 0 ) {
				if ( use_cutpoints_ && cut > 0 ) {
					TR.Debug  <<"selecting loop with cutpoint only" << std::endl;
					found_loops.add_loop( individual_loop[1], individual_loop[individual_loop.size()], cut , 0, false );
				}
				if ( !use_cutpoints_ ) {
					TR.Debug << "not using cutpoint " << std::endl;
					found_loops.add_loop( individual_loop[1], individual_loop[individual_loop.size()], cut , 0, false );
				}
			}
			individual_loop.clear();
			cut = 0;
		}//if statement loop stopped
	}//end for loop through sec strct

	TR <<"loop return: " << found_loops << std::endl;
	protocols::loops::LoopsOP newloops( new protocols::loops::Loops( found_loops ) );
	return newloops;
}


bool
DefineMovableLoops::is_cut( utility::vector1<core::Size> & cut_points, core::Size residue){
	bool res_cut = false;
	for ( core::Size & cut_point : cut_points ) {
		if ( cut_point == residue ) {
			res_cut = true;
		}
	}
	return res_cut;
}

void
DefineMovableLoops::apply( core::pose::Pose & pose ){
	using protocols::loops::Loops;

	if ( secstructure_ == "self" ) {
		secstructure_ = pose.secstruct();
		TR<<"extracting secondary structure from input pose" <<std::endl;
	}

	utility::vector1< core::Size > chains_local = chains_;
	if ( chains_local.empty() ) {
		TR<<"no chains specified, defaulting to use all chains" << std::endl;
		for ( core::Size chain = 1; chain <= pose.conformation().num_chains(); ++chain ) {
			chains_local.push_back( chain );
		}
	}

	core::Size residues = 0;
	//ensure that the residues specified are covered by the secondary structure input
	for ( core::Size chain: chains_local ) {
		residues += pose.split_by_chain( chain )->size();
		TR.Debug <<"residues to compare: "<<residues <<std::endl;
	}
	TR << pose.fold_tree() <<std::endl;
	TR.Debug <<"residues " <<residues <<" ss assignment: "<< secstructure_.length();

	if ( residues != secstructure_.length() ) {
		TR <<"residues vs " <<residues <<" ss assignment: "<< secstructure_.size() << std::endl;
		utility_exit_with_message("input residues under considerations do not agree with the number of secondary strcutres assignments");
	}

	//define offset points, as in which residue to start searching loops
	core::Size start_res = pose.conformation().chain_begin( chains_local[1] );
	//end point, at the last chain
	core::Size stop_res = pose.conformation().chain_end( chains_local[chains_local.size()] );

	//for debugging
	if ( secstructure_.length() != stop_res - start_res + 1 ) {
		TR.Debug << "secstr lenght: " << secstructure_.length() << " stop: " << stop_res << " start: " << start_res << std::endl;
		utility_exit_with_message("secondary structure length does not equal the start and stop of the chain!!!" );
		//allow blueprint incorporation above.... todo
	}

	Loops seeds( parse_seeds(pose, seed_vector_ ) );
	TR.Debug <<"start searching: "<<start_res <<" stop searching: " << stop_res <<std::endl;
	loops_ =  find_loops( pose, secstructure_ , start_res , seeds );///make accessor for seeds
	TR <<"loops " << *loops_ <<std::endl;

}


void
DefineMovableLoops::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	TR<<"DefineMovableLoops has been instantiated"<<std::endl;

	//adding the LoopOP to the data map
	loops_ = utility::pointer::make_shared< protocols::loops::Loops >();
	data.add( "loops", "found_loops", loops_ );

	chainbreakweights_ = tag->getOption< bool >("add_chainbreakterm" , true );

	//get secondary structure either from input template, a string in the xml or through the
	use_cutpoints_ = tag->getOption< bool >( "cutpoint_based" , true );

	if ( tag->hasOption("secstrct") ) {
		secstructure_ = tag->getOption< std::string > ( "secstrct" );
		TR<<"getting secstructure from a string" <<std::endl;
		if ( secstructure_ == "self" ) {
			TR<<"extracting secondary structure from input pose" <<std::endl;
		}
	}

	if ( tag->hasOption( "template_pdb" ) ) {
		std::string const template_pdb_fname( tag->getOption< std::string >( "template_pdb" ));
		template_pdb_ = utility::pointer::make_shared< core::pose::Pose >() ;
		core::import_pose::pose_from_file( *template_pdb_, template_pdb_fname , core::import_pose::PDB_file);
		TR<<"read in a template pdb with " <<template_pdb_->size() <<"residues"<<std::endl;
		core::scoring::dssp::Dssp dssp( *template_pdb_ );
		dssp.insert_ss_into_pose( *template_pdb_ );
		for ( core::Size res = 1 ;  res <= template_pdb_->size(); ++res ) secstructure_ += template_pdb_->secstruct( res );
		secstructure_ = template_pdb_->secstruct();
		TR << secstructure_ << std::endl;
	}

	//command line option overwrites parser input for secondary structure string
	/*
	if (option[ OptionKeys::in::file::psipred_ss2 ].user() ) {
	std::string filename( option[ OptionKeys::in::file::psipred_ss2 ]().name() );
	utility::vector1< char > secstructs;
	utility::io::izstream data( filename );
	//core::pose::util.cc has a method to parse it best
	}*/

	if ( tag->hasOption( "chain_num" ) ) {
		TR<<"NOTE: chains have to be consecutive" << std::endl;
		std::string chain_val( tag->getOption< std::string >( "chain_num" ) );
		utility::vector1< std::string > const chain_keys( utility::string_split( chain_val, ',' ) );
		for ( std::string const & key : chain_keys ) {
			core::Size n;
			std::istringstream ss( key );
			ss >> n;
			chains_.push_back( n );
			TR<<"adding chain "<<key<<std::endl;
		}
	}

	/// read input seeds
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );
	for ( TagCOP btag : branch_tags ) {

		if ( btag->getName() == "Seeds" ) { //need an assertion for the presence of these or at least for the option file
			std::string const beginS( btag->getOption<std::string>( "begin" ) );
			std::string const endS( btag->getOption<std::string>( "end" ) );
			std::pair <std::string,std::string> seedpair;
			seedpair.first  = beginS;
			TR.Debug <<"parsing seeds: " << beginS << " " <<endS <<std::endl;
			seedpair.second = endS;
			seed_vector_.push_back( seedpair );
		}//end seeds
	}//end b-tags
}//end parse_my_tag

std::string DefineMovableLoops::get_name() const {
	return mover_name();
}

std::string DefineMovableLoops::mover_name() {
	return "DefineMovableLoops";
}

void DefineMovableLoops::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default("add_chainbreakterm", xsct_rosetta_bool, "Penalize chain breaks.", "1")
		+ XMLSchemaAttribute::attribute_w_default("cutpoint_based", xsct_rosetta_bool, "Use cutpoint in loop modeling.", "1")
		+ XMLSchemaAttribute("secstrct", xs_string, "File containg sec struct assignments for all residues in \"chain_num\"")
		+ XMLSchemaAttribute("template_pdb", xs_string, "Template pdb file")
		+ XMLSchemaAttribute("chain_num", xs_string, "Comma-separated list of chain IDs. NOTE: chains have to be consecutive");

	AttributeList subelement_attributes;
	subelement_attributes
		+ XMLSchemaAttribute::required_attribute("begin", xs_string, "First residue of seed fragment.")
		+ XMLSchemaAttribute::required_attribute("end", xs_string, "Last residue of seed fragment.");
	XMLSchemaSimpleSubelementList subelement_list;
	subelement_list.add_simple_subelement("Seeds", subelement_attributes, "Defines a seed fragment.");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(),
		"Specific to seeded abinitio protocol. Finds loops after seeded abinitio that should be closed.", attlist, subelement_list );
}

std::string DefineMovableLoopsCreator::keyname() const {
	return DefineMovableLoops::mover_name();
}

protocols::moves::MoverOP
DefineMovableLoopsCreator::create_mover() const {
	return utility::pointer::make_shared< DefineMovableLoops >();
}

void DefineMovableLoopsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DefineMovableLoops::provide_xml_schema( xsd );
}

}
}
