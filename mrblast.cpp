////////////////////////////////////////////////////////////////////////////////
//
//  Parallelizing BLAST on MR-MPI
//
//  Author: Seung-Jin Sul
//          (ssul@jcvi.org)
//
//  Revisions
//
//      09.01.2010  Start!
//
//      09.05.2010  blast call done
//
//      09.07.2010  Adding KV in a mapper and aggregating the vecResults
//
//      09.09.2010  Start v2. Decide to run several MR-MPI jobs and each
//                  job searches on a specific DB chunk.
//
//      09.10.2010  v2 done!
//                  Add print2file member function to KeyValue class.
//
//      09.15.2010  v3 with NCBI C++ Toolkit started.
//
//      09.20.2010  Running version is done. ASN.1 parsing should be done.
//
////////////////////////////////////////////////////////////////////////////////

//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.

/*
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

void    set_default_opts(CRef<CBlastOptionsHandle> optsHandle);
void    run_blast(int itask, char *file, KeyValue *kv, void *ptr);
int     key_compare(char *p1, int len1, char *p2, int len2);
int     key_compare2(char *p1, int len1, char *p2, int len2);

   
 
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

    //if (argc != 4) {
    //if (MPI_myId == 0) printf("Syntax: mpirun -np n mrblast masterFileName chunkName outFileName\n");
    //MPI_Abort(MPI_COMM_WORLD, 1);
    //}

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

    ///
    /// Blast search
    ///
    MapReduce *mr2 = new MapReduce(MPI_COMM_WORLD);
    mr2->verbosity = 0;
    mr2->timer = 1;
    mr2->mapstyle = 2; /// master/slave mode

    int nvecRes;
    char *masterFileName = argv[1];

    giftBox gf;
    gf.myId = MPI_myId;
    gf.numProc = MPI_nProcs;
    gf.pName = MPI_procName;
    gf.dbChunkName = argv[2];

    // Execute main application function
    //int newArgc = 7;
    //char *newArgv[] = {"./mrblast", "-program", "blastn", "-db", "nt.00", "-in", "test2.query"};
    //CBlastDemoApplication demo;
    //demo.AppMain(newArgc, newArgv);
    //CBlastDemoApplication().AppMain(newArgc, newArgv);

    //gf.cbDemo = &demo;
    nvecRes = mr2->map(masterFileName, &run_blast, &gf);
    cout << "### [Node " << MPI_myId << "] nvecRes = " << nvecRes << endl;
    ////mr2->print(-1, 1, 5, 5);

    MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////



    delete mr2;
    MPI::Finalize();

}
 

/* ---------------------------------------------------------------------- */
void run_blast(int itask, char *file, KeyValue *kv, void *ptr)
/* ---------------------------------------------------------------------- */
{
    giftBox *gf = (giftBox *) ptr;
    //CBlastDemoApplication* demo = gf->cbDemo;

    ///
    /// Make args, argv for AppMain
    /// ex) blast -program blastn -db nt.00 -in test.query
    ///
    //int newArgc = 7;
    //char *newArgv[] = {"./mrblast", "-program", "blastn", "-db", "nt.00", "-in", "test2.query"};
    //CBlastDemoApplication().AppMain(newArgc, newArgv);
    //demo->AppMain(newArgc, newArgv);

    ///
    /// Make a option handle and CBlastInputSource
    ///
    
    EProgram program = ProgramNameToEnum("blastn");
    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));
    //CRef<CBlastOptionsHandle> optsHandle(CBlastOptionsFactory::Create(eBlastn));
    opts->Validate();
    
    /// DEBUG
    //opts->GetOptions().DebugDumpText(cerr, "opts", 1);
    
    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
        throw std::runtime_error("Could not initialize object manager");
    }
    
    bool isProtein = false;
    SDataLoaderConfig dlconfig(isProtein);
    CBlastInputSourceConfig iconfig(dlconfig);
    //string qFileName = "test2.query";   
    //ifstream queryFile(qFileName.c_str()); 
    ifstream queryFile(file); 
    //CBlastFastaInputSource fasta_input(qFileName, iconfig);
    
    /// Use query string instead of file input
    string header, seq;
    vector<string> vecSeq;
    unsigned int maxQueryStack = 1;
    
    ///
    /// Target DB setting
    ///
    cout << "dbChunkName = " << gf->dbChunkName << endl;
    string dbChunkName(gf->dbChunkName);
    const CSearchDatabase target_db(dbChunkName, CSearchDatabase::eBlastDbIsNucleotide);
    
    while (!getline(queryFile, header).eof()) {
        getline(queryFile, seq);
        if (header.length() == 0 || seq.length() == 0) {
            fprintf(stderr, "### ERROR: [Node %d]: %s, itask = %d, file = %s, Seq read error.\n",
                    gf->myId, gf->pName, itask, file);
            exit(0);
        }
        fprintf(stdout, "### INFO: [Node %d]: %s, itask = %d, %s\n",
                gf->myId, gf->pName, itask, header.c_str());
    
        string query = header + '\n' + seq;
        
        int j = 1;
        while (j < maxQueryStack && !getline(queryFile, header).eof()) {
            getline(queryFile, seq);
            fprintf(stdout, "### INFO: [Node %d]: %s, itask = %d, %s, stacked = %d\n",
                    gf->myId, gf->pName, itask, header.c_str(), j+1);
            query = query + '\n' + header + '\n' + seq;
            j++;
        }        
        
        ///
        /// Set queries as fasta input
        ///
        CBlastFastaInputSource fasta_input(query, iconfig);
        CBlastInput blastInput(&fasta_input);    
        CScope scope(*objmgr);
        TSeqLocVector queryLoc = blastInput.GetAllSeqLocs(scope);
        CRef<IQueryFactory> queryFactory(new CObjMgr_QueryFactory(queryLoc));
        
        ///
        /// Target DB setting
        ///
        //cout << "dbChunkName = " << gf->dbChunkName << endl;
        //string dbChunkName(gf->dbChunkName);
        //const CSearchDatabase target_db(dbChunkName, CSearchDatabase::eBlastDbIsNucleotide);

        ///
        /// Run blast
        ///    
        CLocalBlast blaster(queryFactory, opts, target_db);
        CSearchResultSet results = *blaster.Run();

        ///
        /// Get the results
        ///

        /// Get warning messages.
        for (unsigned int i = 0; i < results.GetNumResults(); i++) {
            TQueryMessages messages = results[i].GetErrors(eBlastSevWarning);
            if (messages.size() > 0) {
                CConstRef<CSeq_id> seq_id = results[i].GetSeqId();
                if (seq_id.NotEmpty())
                    cerr << "ID: " << seq_id->AsFastaString() << endl;
                else
                    cerr << "ID: " << "Unknown" << endl;

                ITERATE(vector<CRef<CSearchMessage> >, it, messages)
                cerr << (*it)->GetMessage() << endl;
            }
        }

        cout << "results.GetNumResults() = " << results.GetNumResults() << endl;
        for (unsigned int i = 0; i < results.GetNumResults(); i++) {

            CConstRef<CSeq_id> seq_id = results[i].GetSeqId();
            cout << "qID = " << seq_id->AsFastaString() << endl;

            CConstRef<CSeq_align_set> sas = results[i].GetSeqAlign();
            cout << MSerial_AsnText << *sas;
        }
    } // END WHILE

    cout << "Done!\n";
}

/* ---------------------------------------------------------------------- */
void set_default_opts(CRef<CBlastOptionsHandle> optsHandle)
/* ---------------------------------------------------------------------- */
{
    //optsHandle->SetEvalueThreshold(0);
    //optsHandle->SetMatchReward(0);
    //optsHandle->SetMismatchPenalty(0);
    //optsHandle->SetMatrixName("BLOSUM62");
    //optsHandle->SetCutoffScore();
    if (CBlastNucleotideOptionsHandle* nuclHandle =
        dynamic_cast<CBlastNucleotideOptionsHandle*>(&*optsHandle)) {

        //nuclHandle->SetMatchReward(0);
        //nuclHandle->SetMismatchPenalty(0);
        nuclHandle->SetMatrixName("BLOSUM62");
        //nuclHandle->SetCutoffScore();
    }

    return;
}

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

