// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief declaration of implementation class for abstract class Residue
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_AtomICoor_hh
#define INCLUDED_core_chemical_AtomICoor_hh


// Unit headers
#include <core/chemical/AtomICoor.fwd.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

#include <core/id/AtomID.fwd.hh>


// C++ headers


namespace core {
namespace chemical {

/// @brief Atom 's ID in internal coordinates in a ResidueType
class ICoorAtomID {
public:
	typedef conformation::Residue Residue;
	typedef conformation::Conformation Conformation;

public:


public:
	/// @brief default constructor
	ICoorAtomID();

	/// @brief construct ICoorAtomID by atom name and its ResidueType
	ICoorAtomID(
		std::string name,
		ResidueType const & rsd_type
	);

public:
	/// @brief get ICoorAtomID atomno
	Size
	atomno() const
	{
		return atomno_;
	}

	/// @brief set ICoorAtomID atomno
	void
	atomno( int const atomno_in )
	{
		atomno_ = atomno_in;
	}

	/// @brief get ICoordAtomID type
	ICoordAtomIDType const &
	type() const
	{
		return type_;
	}


	bool
	is_internal() const
	{
		return ( type_ == ICoordAtomIDType::INTERNAL );
	}


	bool
	is_polymer_lower() const
	{
		return ( type_ == ICoordAtomIDType::POLYMER_LOWER );
	}


	bool
	is_polymer_upper() const
	{
		return ( type_ == ICoordAtomIDType::POLYMER_UPPER );
	}

	/// @brief Returns true if this is the specified connection id
	///
	bool
	is_connect( Size const connid ) const
	{
		return ( type_ == ICoordAtomIDType::CONNECT && atomno_ == connid );
	}

	/// @brief Returns true if this is a connection.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool
	is_connect( ) const
	{
		return ( type_ == ICoordAtomIDType::CONNECT );
	}

	/// @brief Returns the string representation which will build this ICoorAtomID
	/// (e.g. atom name, UPPER, LOWER, CONN*
	std::string
	name( ResidueType const & rt ) const;

public:


	/// @brief What is the coordinates corresponding to this ICoorAtomID,
	/// for the given residue and conformation
	Vector const &
	xyz( Residue const & rsd, Conformation const & conformation ) const;


	/// @brief What is the coordinates corresponding to this ICoorAtomID,
	/// for the given idealized ResidueType
	Vector
	xyz( ResidueType const & rsd_type ) const;

	/// @brief WARNING: Slightly dangerous function intended for black magic use only.
	///    Rebuilds atom location from stub atoms. If stub atom are not internal atoms, their
	///    location will be rebuilt from their residue stub atom's locations, as opposed to being
	///    retrieved from connected residues via a conformation.
	Vector
	xyz( conformation::Residue const & rsd ) const;

	/// @brief This ICoorAtomID (for the given residue) corresponds to which id::AtomID in the conformation?
	id::AtomID
	atom_id( Residue const & rsd, Conformation const & conformation ) const;

	/// @brief Can valid coordinates be built for this ICoorAtomID,
	/// given the residue and conformation?
	bool
	buildable( Residue const & rsd, Conformation const & conformation ) const;

	bool operator==( ICoorAtomID const & other ) const;
	bool operator!=( ICoorAtomID const & other ) const;

private:

	/// atom's "connection" type
	ICoordAtomIDType type_;
	/// atom's index number
	Size atomno_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief A basic class containing info of internal coordinates needed for building an atom within a ResidueType
/**
In atom tree, each atom is defined by its internal coordinates, which include a bond distance,
a bond angle and a torsion angle. Of course, all these internal coordinates are only meaningful
in the context of three reference (stub) atoms. AtomICoor information is stored in the residue param
files and some terms are defined as following:\n
- bond distance d_ is that between the atom to be built (child) and stub_atom1 (parent)
- bond angle theta_ is that defined by child-parent-stub2(angle)
- torsion angle phi_ is that defined by child-parent-stub2-stub3(torsion)
*/
class AtomICoor {
public:
	/// @brief default constructor
	AtomICoor();
	/// @brief constructor
	AtomICoor(
		std::string const & built_atom_name,
		Real const phi_in,
		Real const theta_in,
		Real const d_in,
		std::string const & stub_atom1_name,
		std::string const & stub_atom2_name,
		std::string const & stub_atom3_name,
		ResidueType const & rsd_type
	);

	AtomICoor(
		std::string const & built_atom_name,
		Real const phi_in,
		Real const theta_in,
		Real const d_in,
		ICoorAtomID const & stub_atom1,
		ICoorAtomID const & stub_atom2,
		ICoorAtomID const & stub_atom3,
		ResidueType const & rsd_type
	);

public:
	/// @brief accessor to stub_atom1 ICoorAtomID
	Real
	phi() const
	{
		return phi_;
	}


	Real
	theta() const
	{
		return theta_;
	}


	Real
	d() const
	{
		return d_;
	}


	ICoorAtomID const &
	stub_atom1() const
	{
		return stub_atom1_;
	}

	/// @brief accessor to stub_atom2 ICoorAtomID
	ICoorAtomID const &
	stub_atom2() const
	{
		return stub_atom2_;
	}

	/// accessor to stub_atom3 ICoorAtomID
	ICoorAtomID const &
	stub_atom3() const
	{
		return stub_atom3_;
	}


	bool
	is_internal() const
	{
		return ( stub_atom1_.is_internal() && stub_atom2_.is_internal() && stub_atom3_.is_internal() );
	}

	/// @brief Given a Boolean vector corresponding to atom indices, determine whether any of the stub
	/// indices depends on one of the atoms that are "true".
	/// @details Called by ResidueType::update_polymer_dependent_groups() when figuring out which atoms depend on atoms that
	/// depend on a connection.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool
	depends_on_a_true_index( utility::vector1< bool > const &atomvect ) const {
		for ( core::Size i(1), imax(atomvect.size()); i<=imax; ++i ) {
			if ( atomvect[i] ) {
				if ( stub_atom1_.atomno() == i || stub_atom2_.atomno() == i || stub_atom3_.atomno() == i ) return true;
			}
		}
		return false;
	}

	/// @brief Is this icoor immediately dependent on LOWER?
	/// @details Returns true if stub atom 1, 2, or 3 is LOWER.
	bool
	depends_on_polymer_lower() const
	{
		return ( stub_atom1_.is_polymer_lower() || stub_atom2_.is_polymer_lower() || stub_atom3_.is_polymer_lower() );
	}


	/// @brief Is this icoor immediately dependent on UPPER?
	/// @details Returns true if stub atom 1, 2, or 3 is UPPER.
	bool
	depends_on_polymer_upper() const
	{
		return ( stub_atom1_.is_polymer_upper() || stub_atom2_.is_polymer_upper() || stub_atom3_.is_polymer_upper() );
	}

	/// @brief Returns true if any of the stub atoms is the specified connection ID.
	///
	bool
	depends_on_residue_connection( Size const connid ) const
	{
		return ( stub_atom1_.is_connect( connid ) ||
			stub_atom2_.is_connect( connid ) ||
			stub_atom3_.is_connect( connid ) );
	}

	/// @brief Returns true if any of the stub atoms is a connection ID.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool
	depends_on_residue_connection( ) const
	{
		return ( stub_atom1_.is_connect( ) ||
			stub_atom2_.is_connect( ) ||
			stub_atom3_.is_connect( ) );
	}


	/// @brief accessor to stub_atom ICoorAtomID
	ICoorAtomID &
	stub_atom( int const atm )
	{
		switch( atm ) {
		case 1 : return stub_atom1_;
		case 2 : return stub_atom2_;
		case 3 : return stub_atom3_;
		}
		utility_exit_with_message( "ICoorAtomID::stub_atom should be 1--3" );
		return stub_atom1_;
	}

	/// @brief constant accessor to stub_atom ICoorAtomID
	ICoorAtomID const &
	stub_atom( int const atm ) const
	{
		switch( atm ) {
		case 1 : return stub_atom1_;
		case 2 : return stub_atom2_;
		case 3 : return stub_atom3_;
		}
		utility_exit_with_message( "ICoorAtomID::stub_atom should be 1--3" );
		return stub_atom1_;
	}

	/// @brief The name of the atom being built by this icoor
	std::string const & built_atom() const
	{
		return built_atom_;
	}


	void show( std::ostream & ) const;

public:

	Vector
	build(
		conformation::Residue const & rsd,
		conformation::Conformation const & conformation
	) const;

	Vector
	build(
		ResidueType const & rsd_type
	) const;

	/// @brief WARNING: Slightly dangerous function intended for black magic use only.
	///    Rebuilds atom location from stub atoms. If stub atom are not internal atoms, their
	///    location will be rebuilt from their residue stub atom's locations, as opposed to being
	///    retrieved from connected residues via a conformation.
	Vector
	build( conformation::Residue const & rsd ) const;

private:

	/// @brief Internal, centralized implementation of the build() functions.
	Vector
	build(
		Vector v1,
		Vector v2,
		Vector v3
	) const;

private:

	std::string built_atom_;
	Real phi_;
	Real theta_;
	Real d_;
	ICoorAtomID stub_atom1_;
	ICoorAtomID stub_atom2_;
	ICoorAtomID stub_atom3_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief print the icoord table starting at the root.
void pretty_print_atomicoor(std::ostream & out, ResidueType const & rsd_type);

/// @brief print the icoord table starting with the given icoord record
void pretty_print_atomicoor(std::ostream & out, AtomICoor const & start, ResidueType const & rsd_type, core::Size indent=0);

std::ostream & operator <<( std::ostream & out, ICoordAtomIDType type );

} // chemical
} // core


#endif
