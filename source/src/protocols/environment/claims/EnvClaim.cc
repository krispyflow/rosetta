// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/environment/claims/EnvClaim.cc
/// @brief  Abstract class for implementing an EnvClaim (e.g. TorsionClaim).
/// @author Justin R. Porter

// Unit Headers
#include <protocols/environment/claims/EnvClaim.hh>

// Package Headers
#include <protocols/environment/ClientMover.hh>
#include <core/environment/SequenceAnnotation.hh>


#include <protocols/environment/claims/CutBiasClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>
#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/VirtResClaim.hh>
#include <protocols/environment/claims/XYZClaim.hh>

// Project Headers

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

#include <utility>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

#include <boost/algorithm/string.hpp>

//// C++ headers
#include <iostream>

#include <core/select/residue_selector/ResidueSelector.hh> // AUTO IWYU For ResidueSelector

// option key includes


static basic::Tracer tr( "protocols.environment.EnvClaim", basic::t_info );

namespace protocols {
namespace environment {
namespace claims {

EnvClaimOP EnvClaim::make_claim( std::string const& name,
	ClientMoverOP owner,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap ) {
	if      ( name == "CutBiasClaim" ) return utility::pointer::make_shared< CutBiasClaim >( owner, tag, datamap );
	else if ( name == "JumpClaim" )    return utility::pointer::make_shared< JumpClaim >( owner, tag, datamap );
	else if ( name == "TorsionClaim" ) return utility::pointer::make_shared< TorsionClaim >( owner, tag, datamap );
	else if ( name == "VirtResClaim" )  return utility::pointer::make_shared< VirtResClaim >( owner, tag, datamap );
	else if ( name == "XYZClaim" )     return utility::pointer::make_shared< XYZClaim >( owner, tag, datamap );
	else {
		tr << "NOTE: The VrtResClaim is now called VirtResClaim. Please alter your scripts accordingly." << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "'" + name + "' is not a known EnvClaim type." );
	}
}

std::string
EnvClaim::envclaim_ct_namer( std::string tag_name ){
	return "envclaim_" + tag_name + "_complex_type";
}

std::string
EnvClaim::envclaim_group_name(){
	return "envclaim";
}


void
EnvClaim::define_envclaim_schema_group( utility::tag::XMLSchemaDefinition & xsd ){
	//NOTE: Since there is no factory system or creators for EnvClaims, these are just defined manually here.
	//Ideally these should be set up with a factory system and creators registered in protocols/init
	using namespace utility::tag;
	CutBiasClaim::provide_xml_schema( xsd );
	JumpClaim::provide_xml_schema( xsd );
	TorsionClaim::provide_xml_schema( xsd );
	VirtResClaim::provide_xml_schema( xsd );
	XYZClaim::provide_xml_schema( xsd );


	XMLSchemaElementOP cutbias_element( new XMLSchemaElement );
	cutbias_element->name( CutBiasClaim::class_name() );
	cutbias_element->type_name( envclaim_ct_namer( CutBiasClaim::class_name() ) );

	XMLSchemaElementOP jump_element( new XMLSchemaElement );
	jump_element->name( JumpClaim::class_name() );
	jump_element->type_name( envclaim_ct_namer( JumpClaim::class_name() ) );

	XMLSchemaElementOP torsion_element( new XMLSchemaElement );
	torsion_element->name( TorsionClaim::class_name() );
	torsion_element->type_name( envclaim_ct_namer( TorsionClaim::class_name() ) );

	XMLSchemaElementOP vrt_element( new XMLSchemaElement );
	vrt_element->name( VirtResClaim::class_name() );
	vrt_element->type_name( envclaim_ct_namer( VirtResClaim::class_name() ) );

	XMLSchemaElementOP xyz_element( new XMLSchemaElement );
	xyz_element->name( XYZClaim::class_name() );
	xyz_element->type_name( envclaim_ct_namer( XYZClaim::class_name() ) );

	XMLSchemaModelGroupOP envclaim_choice( new XMLSchemaModelGroup );
	envclaim_choice->type( xsmgt_choice );
	envclaim_choice->append_particle( cutbias_element );
	envclaim_choice->append_particle( jump_element );
	envclaim_choice->append_particle( torsion_element );
	envclaim_choice->append_particle( vrt_element );
	envclaim_choice->append_particle( xyz_element );

	XMLSchemaModelGroup envclaim_group;
	envclaim_group.group_name( envclaim_group_name() );
	envclaim_group.append_particle( envclaim_choice );
	xsd.add_top_level_element( envclaim_group );
}











bool EnvClaim::is_claim( std::string const& name ) {
	if      ( name == "CutBiasClaim" ) return true;
	else if ( name == "JumpClaim" )    return true;
	else if ( name == "TorsionClaim" ) return true;
	else if ( name == "VirtResClaim" )  return true;
	else if ( name == "XYZClaim" )     return true;
	else return false;
}

EnvClaim::EnvClaim( ClientMoverOP owner ):
	VirtualBase(),
	claim_source_(std::move( owner ))
{}

/// @details Auto-generated virtual destructor
EnvClaim::~EnvClaim() = default;

void EnvClaim::show( std::ostream& os ) const {
	os << "owned by, " << owner()->type() << ";";
}

DOFElement EnvClaim::wrap_dof_id( core::id::DOF_ID const& id ) const {
	DOFElement e;
	e.id = id;

	return e;
}

ClientMoverOP EnvClaim::owner() const {
	return claim_source_;
}

ControlStrength EnvClaim::parse_ctrl_str( std::string const& str ) const {
	std::string lower = str;
	boost::algorithm::to_lower( lower );

	if ( lower == "does_not_control" ) {
		return DOES_NOT_CONTROL;
	} else if ( lower == "can_control" ) {
		return CAN_CONTROL;
	} else if ( lower == "must_control" ) {
		return MUST_CONTROL;
	} else if ( lower == "exclusive" ) {
		return EXCLUSIVE;
	} else {
		throw CREATE_EXCEPTION(utility::excn::BadInput,  "The initialization strength '" + str +
			"' is not recognized." );
	}
}

void EnvClaim::annotate( core::pose::Pose const& pose, core::environment::SequenceAnnotationOP ann ) const {
	for ( AnnotatingSelectors::value_type const & pair : selector_list_ ) {
		std::string const& label = pair.first;
		ResidueSelectorCOP selector = pair.second;

		utility::vector1< bool > subset = selector->apply( pose );

		utility::vector1< core::Size > trues;
		for ( core::Size i = 1; i <= subset.size(); ++i ) {
			if ( subset[i] ) {
				trues.push_back( i );
			}
		}

		try{
			ann->add_seq_label( label, trues );
		} catch ( utility::excn::KeyError& e ) {
			std::ostringstream ss;
			ss << "While " << *this << " was annotating the pose for broking, the " << selector->get_name()
				<< "Selector produced a conflicting residue selection.";
			e.add_msg( ss.str() );
			throw;// e;
		}
	}
}

void EnvClaim::queue_for_annotation( std::string const& label, ResidueSelectorCOP selector ) {
	if ( selector_list_.find( label ) == selector_list_.end() ||            // no selector with label 'label' exists
			selector_list_.find( label )->second.get() == selector.get() ) {    // selector with label 'label' *is* 'selector'
		selector_list_[ label ] = selector;
	} else {
		std::ostringstream ss;
		ss << "In claim " << *this << ", the label '" << label << "' was to be used for the selection of a "
			<< selector->get_name() << "Selector, but that label was already used by an existing "
			<< selector_list_[label]->get_name() << "Selector.";
		if ( selector->get_name() == selector_list_[label]->get_name() ) {
			ss << " They aren't the same ResidueSelector.";  // Make sure error is clear that selector-sameness has be checked already.
		}
		throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
	}
}

extern std::ostream& operator<<( std::ostream& os, EnvClaim const& claim ) {
	claim.show( os );
	return os;
}

extern std::ostream& operator<<( std::ostream& os, EnvClaims const& claims ) {
	for ( auto const & claim : claims ) {
		if ( claim ) {
			os << *claim << "\n";
		} else {
			os << "No-Claim\n";
		}
	}
	return os;
}

} //claims
} //environment
} //protocols
