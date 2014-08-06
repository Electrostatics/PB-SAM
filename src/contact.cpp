#include "contact.h"

REAL CContact::SEPDIST = 2.0;

/******************************************************************/
/******************************************************************/
/**
*  Constructor to create MolContact class
******************************************************************/
CMolContact::CMolContact( char* fname ) : m_bDocked(false)
{
	m_contactlist.clear(  );
	
	cout <<"Reading molcontact file "<<fname<<endl;
	cout <<"Using sepdist criteria = "<<CContact::SEPDIST<<endl;
	
	// File handling
	ifstream fin( fname );
	if ( !fin.is_open( ))
	{
		// Unable to open file, die
		cout << "Could not open file " << fname << endl;
		exit( 0 );
	}
	
	fin >> m_mol2type;
	
	fin >> m_ncontact;
	cout << m_mol2type<<" "<<m_ncontact<<endl;
	
	while ( !fin.eof( ))
	{
		int k1,k2;
		fin >>k1>>k2; 
		m_contactlist.push_back(  CContact(k1, k2 ));
	}
	
	fin.close(  );
} // end CMolContact : m_bDocked( false )

/******************************************************************/
/******************************************************************/
/**
*  Constructor to create MolContact class
******************************************************************/
CMolContact::CMolContact( int mol2type, int ncontact, const vector<CContact> &clist ) 
: m_mol2type( mol2type ), m_ncontact(ncontact), m_contactlist(clist), m_bDocked(false)
{}

/******************************************************************/
/******************************************************************/
/**
*  read in molcontacts
******************************************************************/
void
CMolContact::readMolContact( const char* fname, int &mol1type, 
														vector<CMolContact> &molcontactlist, double sepdist )
{
	cout <<"Reading molcontact file "<<fname<<endl;
	cout <<"Using sepdist criteria = "<<sepdist<<endl;
	
	ifstream fin( fname );
	if ( !fin.is_open( ))
	{
		cout << "Could not open file " << fname << endl;
		exit( 0 );
	}
	
	int mol2type, ndef;
	fin >> mol1type; 
	fin >> ndef;
	int npair, ncontact;
	int line = 2;
	int tnpair=0;
	while ( true )
	{      
		fin >> mol2type; // what the partner molecule should be
		if( fin.eof( )) break;
		
		fin >> ncontact; // number of contacts for docked criteria
		fin >> npair; // number of contact pairs
		line += 3;
		tnpair += npair;
		
		vector<CContact> clist;
		for( int n=0; n<npair;n++ )
		{
			int k1,k2;
			
			fin >>k1>>k2; 
			clist.push_back( CContact(k1, k2, sepdist ));
		}
		
		if( clist.size( ) == npair) 
		{
			line += npair;
			molcontactlist.push_back(  CMolContact(mol2type,ncontact,clist ) );
		}
	}
	
	cout <<"mol definitions read "<<molcontactlist.size(  )<<endl;
	fin.close(  );
}

/******************************************************************/
/******************************************************************/
/**
*  read molcontacts from file
******************************************************************/
/*void
CMolContact::readMolContactWithDist( const char* fname, int &mol1type, 
																		vector<CMolContact> &molcontactlist )
{
	cout <<"Reading molcontact file "<<fname<<endl;
	
	ifstream fin( fname );
	if ( !fin.is_open( ))
	{
		cout << "Could not open file " << fname << endl;
		exit( 0 );
	}
	
	int mol2type, ndef;
	fin >> mol1type; 
	fin >> ndef; // i.e. number of unique binding partners
	int npair, ncontact;
	int line = 2;
	int tnpair=0;
	
	while ( !fin.eof( ))
	{      
		fin >> mol2type; // what the partner molecule should be
		if( fin.eof( )) break;
		fin >> ncontact; // number of contacts for docked criteria
		fin >> npair; // number of contact pairs
		line += 3;
		tnpair += npair;
		
		vector<CContact> clist;
		for( int n=0; n<npair;n++ )
		{
			int k1,k2;
			double dist;
			fin >>k1>>k2>>dist; 
			clist.push_back( CContact(k1, k2, dist ));
		}
		
		if( clist.size( ) == npair) 
		{
			line += npair;
			molcontactlist.push_back(  CMolContact(mol2type,ncontact,clist ) );
		}
	}
	
	cout <<"mol definitions read "<<molcontactlist.size(  )<<endl;
	fin.close(  );
}
*/

/******************************************************************/
/******************************************************************/
/**
*  general, for proteins to bind multiple partners 
******************************************************************/
bool
CMolContact::IsDocked( CMolecule* mol1, CMolecule* mol2, 
											vector<vector<CMolContact> > MOLCONTACTLIST ) 
{
	int m1 = mol1->getMolType(  );
	
	int nDockDef = MOLCONTACTLIST[m1].size(  );
	for( int h=0; h<nDockDef; h++ )
	{
		CMolContact mcon = MOLCONTACTLIST[m1][h];
		
		if( mcon.getMol2Type( )!= mol2->getMolType() )  continue;
		
		int ncon = 0;
		vector<CContact> clist = mcon.getContactList(  );
		for( int k=0; k<clist.size( ) && ncon < mcon.getNContact(); k++)
		{
			if(  CMolecule::getSepDist(mol1, clist[k].getID1( ), mol2, clist[k].getID2())
				 <= clist[k].getDist(  ) ) ncon++;
			
			if(  ncon == mcon.getNContact( )) return true;
		}
	} // end all definitions                                                                         
	
	return false;
}

