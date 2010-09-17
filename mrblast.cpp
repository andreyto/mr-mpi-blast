/*  $Id: blast_demo.cpp 144153 2008-10-27 19:49:01Z vakatov $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Tom Madden
 *
 * File Description:
 *   Sample application for the running a blast search.
 *
 */

/// MPI AND MR-MPI
#include "mpi.h"
#include "./mrmpi/mapreduce.h"
#include "./mrmpi/keyvalue.h"

using namespace MAPREDUCE_NS;
using namespace std;
 
const size_t MAX_STR = 256;


/// NCBI 
#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbiargs.hpp>

#include <objmgr/object_manager.hpp>

#include <objects/seqalign/Seq_align_set.hpp>

#include <algo/blast/api/sseqloc.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/blast_types.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>
#include <algo/blast/api/blast_prot_options.hpp>

#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>

USING_NCBI_SCOPE;
USING_SCOPE(blast);

/// TO COLLECT PROC NAMES ON WHICH EACH TASK IS ASIGNED
struct procName {
    int     myId;
    char*   pName;
};

/// FOR COMMUNICATING WITH MR-MPI
struct giftBox {
    int             myId;
    int             numProc;
    int             numDbChunks;
    int             totalDbChunks;
    char*           pName;
    char*           queryFileName;
    char*           dbChunkName;
    unsigned long   numQuery;
    //multimap<string, unsigned int> *mmapPName;
    vector<string>* vPName;
};

/// MR-MPI CALLS
void    get_node_name(int, KeyValue *, void *);
void    collect_node_names(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void    save_node_names(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);
//void    run_blast(int itask, KeyValue *kv, void *ptr);
void    run_blast2(int itask, char *file, KeyValue *kv, void *ptr);
int     key_compare(char *p1, int len1, char *p2, int len2);
int     key_compare2(char *p1, int len1, char *p2, int len2); 
 

/////////////////////////////////////////////////////////////////////////////
//  CBlastDemoApplication::


class CBlastDemoApplication : public CNcbiApplication
{
private:
    virtual void Init(void);
    virtual int  Run(void);
    virtual void Exit(void);

    void ProcessCommandLineArgs(CRef<CBlastOptionsHandle> opts_handle);

};


/////////////////////////////////////////////////////////////////////////////
//  Init test for all different types of arguments


void CBlastDemoApplication::Init(void)
{
    // Create command-line argument descriptions class
    auto_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);

    // Specify USAGE context
    arg_desc->SetUsageContext(GetArguments().GetProgramBasename(), "BLAST demo program");

    arg_desc->AddKey
        ("program", "ProgramName",
         "One of blastn, megablast, disc_megablast, blastp, blastx, tblastn, tblastx, rpsblast",
         CArgDescriptions::eString);
    arg_desc->SetConstraint
        ("program", &(*new CArgAllow_Strings,
                "blastn", "megablast", "disc_megablast", "blastp", "blastx", "tblastn", "tblastx", "rpsblast"));

    arg_desc->AddDefaultKey
        ("db", "DataBase",
         "This is the name of the database",
         CArgDescriptions::eString, "nr");

    arg_desc->AddDefaultKey("in", "Queryfile",
                        "A file with the query", CArgDescriptions::eInputFile, "stdin");

    arg_desc->AddDefaultKey("out", "Outputfile",
                        "The output file", CArgDescriptions::eOutputFile, "stdout");

    arg_desc->AddDefaultKey("evalue", "evalue",
                        "E-value threshold for saving hits", CArgDescriptions::eDouble, "0");

    arg_desc->AddDefaultKey("penalty", "penalty", "Penalty score for a mismatch",
                            CArgDescriptions::eInteger, "0");

    arg_desc->AddDefaultKey("reward", "reward", "Reward score for a match",
                            CArgDescriptions::eInteger, "0");

    arg_desc->AddDefaultKey("matrix", "matrix", "Scoring matrix name",
                            CArgDescriptions::eString, "BLOSUM62");

    // Setup arg.descriptions for this application
    SetupArgDescriptions(arg_desc.release());
}

/// Modify BLAST options from defaults based upon command-line args.
///
/// @param opts_handle already created CBlastOptionsHandle to modify [in]
void CBlastDemoApplication::ProcessCommandLineArgs(CRef<CBlastOptionsHandle> opts_handle)

{
    CArgs args = GetArgs();

        // Expect value is a supported option for all flavors of BLAST.
        if(args["evalue"].AsDouble())
          opts_handle->SetEvalueThreshold(args["evalue"].AsDouble());
        
        // The first branch is used if the program is blastn or a flavor of megablast
        // as reward and penalty is a valid option.
        //
        // The second branch is used for all other programs except rpsblast as matrix
        // is a valid option for blastp and other programs that perform protein-protein
        // comparisons.
        //
        if (CBlastNucleotideOptionsHandle* nucl_handle =
              dynamic_cast<CBlastNucleotideOptionsHandle*>(&*opts_handle)) {

              if (args["reward"].AsInteger())
                nucl_handle->SetMatchReward(args["reward"].AsInteger());
            
              if (args["penalty"].AsInteger())
                nucl_handle->SetMismatchPenalty(args["penalty"].AsInteger());
        }
        else if (CBlastProteinOptionsHandle* prot_handle =
               dynamic_cast<CBlastProteinOptionsHandle*>(&*opts_handle)) {
              if (args["matrix"]) 
                prot_handle->SetMatrixName(args["matrix"].AsString().c_str());
        }

        return;
}


/////////////////////////////////////////////////////////////////////////////
//  Run test (printout arguments obtained from command-line)


int CBlastDemoApplication::Run(void)
{    
    // Get arguments
    const CArgs& args = GetArgs();

    EProgram program = ProgramNameToEnum(args["program"].AsString());

    bool db_is_aa = (program == eBlastp || program == eBlastx ||
                     program == eRPSBlast || program == eRPSTblastn);

    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));

    ProcessCommandLineArgs(opts);

    opts->Validate();  // Can throw CBlastException::eInvalidOptions for invalid option.


    // This will dump the options to stderr.
    // opts->GetOptions().DebugDumpText(cerr, "opts", 1);

    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
         throw std::runtime_error("Could not initialize object manager");
    }

    const bool is_protein =
        !!Blast_QueryIsProtein(opts->GetOptions().GetProgramType());
    SDataLoaderConfig dlconfig(is_protein);
    CBlastInputSourceConfig iconfig(dlconfig);
    CBlastFastaInputSource fasta_input(args["in"].AsInputFile(), iconfig);
    CScope scope(*objmgr);

    CBlastInput blast_input(&fasta_input);

    TSeqLocVector query_loc = blast_input.GetAllSeqLocs(scope);

    CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc));

    const CSearchDatabase target_db(args["db"].AsString(),
        db_is_aa ? CSearchDatabase::eBlastDbIsProtein : CSearchDatabase::eBlastDbIsNucleotide);

    CLocalBlast blaster(query_factory, opts, target_db);

    CSearchResultSet results = *blaster.Run();

    // Get warning messages.
    for (unsigned int i = 0; i < results.GetNumResults(); i++) 
    {
        TQueryMessages messages = results[i].GetErrors(eBlastSevWarning);
        if (messages.size() > 0)
        {
            CConstRef<CSeq_id> seq_id = results[i].GetSeqId();
            if (seq_id.NotEmpty())
                cerr << "ID: " << seq_id->AsFastaString() << endl;
            else
                cerr << "ID: " << "Unknown" << endl;

            ITERATE(vector<CRef<CSearchMessage> >, it, messages)
                cerr << (*it)->GetMessage() << endl;
        }
    }
    
    CNcbiOstream& out = args["out"].AsOutputFile();

    for (unsigned int i = 0; i < results.GetNumResults(); i++) {
         CConstRef<CSeq_align_set> sas = results[i].GetSeqAlign();
         //out << MSerial_AsnText << *sas;
         out << MSerial_Xml << *sas;
    }
    
    cout << "Done!!\n";

    return 0;
}


/////////////////////////////////////////////////////////////////////////////
//  Cleanup


void CBlastDemoApplication::Exit(void)
{
    SetDiagStream(0);
}


/////////////////////////////////////////////////////////////////////////////
//  MAIN


#ifndef SKIP_DOXYGEN_PROCESSING
int main(int argc, char* argv[])
{  
    /// MPI
    MPI_Init(&argc, &argv);

    char MPI_procName[100];

    int MPI_myId, MPI_nProcs, MPI_length;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myId);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_nProcs);
    MPI_Get_processor_name(MPI_procName, &MPI_length);
    fprintf(stdout, "### [Node %d]: %s \n", MPI_myId, MPI_procName);
    
    ///MR-MPI
    MapReduce *mr = new MapReduce(MPI_COMM_WORLD);
    mr->verbosity = 0;
    mr->timer = 1;    
    
    MPI_Barrier(MPI_COMM_WORLD);
    //double tstart = MPI_Wtime();
    
    ///
    /// GET MPI PROC NAMES
    /// TO ASSIGN A DB CHUNK WITH A NODE SO THAT THE TASKS ON THE NODE
    /// SEARCHES QUERIES ON THE DB CHUNK
    ///
    procName nn;
    nn.myId = MPI_myId;
    nn.pName = MPI_procName;
    int sizeKV = mr->map(MPI_nProcs, &get_node_name, &nn);
    //mr->print(-1, 1, 5, 5);

    /// DEBUG
    //mr->collate(NULL);
    //mr->print(-1, 1, 5, 5);
    //int numUniqKV = mr->reduce(&collect_node_names, (void *)NULL);
    //mr->print(-1, 1, 5, 5);

    MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////

    mr->gather(1); /// KVs FROM TASKs -> A KV IN THE ROOT NODE
    mr->print(-1, 1, 5, 5);

    mr->broadcast(0);
    //mr->print(-1, 1, 5, 5);

    //multimap<string, unsigned int> mmapProcNames;
    vector<string> vecProcName;
    mr->map(mr, &save_node_names, &vecProcName);
    //cout << vecProcName.size() << endl;

    MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////
    
    /// DEBUG
    //if (MPI_myId == 0) {
        ////pair<multimap<string, unsigned int>::iterator, multimap<string, unsigned int>::iterator> ppp;
        ////ppp = mmapProcNames.equal_range("ssul-lx");
        ////for (multimap<string, unsigned int>::iterator it2 = ppp.first; it2 != ppp.second; ++it2) {
            ////cout << "  [" << (*it2).first << ", " << (*it2).second << "]" << endl;
        ////}
        //for (int i = 0; i < vecProcName.size(); ++i)
            //cout << "### Node " << MPI_myId << " name = " << vecProcName[i] << endl;
    //}
    
    delete mr;
    
    
    
    
    
    
    
    MPI::Finalize();
    
    
    
    
    
    // Execute main application function
    return CBlastDemoApplication().AppMain(argc, argv);    
}
#endif /* SKIP_DOXYGEN_PROCESSING */

/* ---------------------------------------------------------------------- */
void get_node_name(int itask, KeyValue *kv, void *ptr)
/* ---------------------------------------------------------------------- */
{
    procName *nn = (procName *) ptr;
    char id[MAX_STR];
    sprintf(id, "%d", nn->myId);
    kv->add(nn->pName, strlen(nn->pName) + 1, id, strlen(id) + 1);
}

/* ---------------------------------------------------------------------- */
//void collect_node_names(char *key, int keybytes, char *multivalue,
//                        int nvalues, int *valuebytes, KeyValue *kv,
//                        void *ptr)
/* ---------------------------------------------------------------------- */
//{
//kv->add(key, keybytes, multivalue, sizeof(multivalue));
//}

/* ---------------------------------------------------------------------- */
void save_node_names(uint64_t itask, char *key, int keybytes, char *value,
                     int valuebytes, KeyValue *kv, void *ptr)
/* ---------------------------------------------------------------------- */
{
    //multimap<string, unsigned int> *nNames = (multimap<string, unsigned int> *) ptr;
    vector<string>* vecNames = (vector<string> *) ptr;
    string pName(key);
    //unsigned int myid = atoi(value);
    //nNames->insert(pair<string, unsigned int>(pName, myid));
    vector<string>::const_iterator loc = find(vecNames->begin(), vecNames->end(), pName);
    if (loc == vecNames->end())
        vecNames->push_back(pName);
}


/// EOF

