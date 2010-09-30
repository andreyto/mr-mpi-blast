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
//      09.21.2010  XML parsing tried but returned back to use raw blast 
//                  output for creating KV
//
//      09.22.2010  Using raw blast output, KV is constructed. Running 
//                  version is done. Start testing on Ranger.
//
//      09.29.2010  Using work items, v4 is done.
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

/// ASN
#include <sstream>

/// blastn importing
//#include <algo/blast/blastinput/blastn_args.hpp>
//#include <algo/blast/format/blast_format.hpp>
//#include "blast_app_util.hpp"

/// Tabular
//#include <algo/blast/format/blast_format.hpp>
//#include <algo/blast/format/blastxml_format.hpp>
//#include <algo/blast/format/data4xmlformat.hpp>       /* NCBI_FAKE_WARNING */
//#include <algo/blast/format/build_archive.hpp>
//#include <objects/seq/Seq_annot.hpp>
//#include <objects/general/User_object.hpp>
//#include <objects/general/User_field.hpp>
//#include <algo/blast/core/blast_stat.h>
//#include <corelib/ncbiutil.hpp>                 // for FindBestChoice

/// XML
//#include <misc/xmlwrapp/xmlwrapp.hpp>           // for XML::
//#include <misc/xmlwrapp/xsltwrapp.hpp>

/// CSeq_align for output
#include <objects/seqalign/Seq_align.hpp>
#include <util/range.hpp>


USING_NCBI_SCOPE;
USING_SCOPE(blast);
//USING_SCOPE(objects);
//USING_SCOPE(align_format);
//USING_SCOPE(xml);

const unsigned int MAXSTR         = 256;
const unsigned int MAXQUERYNUM    = 1000;

const unsigned int QUERY          = 0;
const unsigned int SUBJECT        = 1;

/// TO COLLECT PROC NAMES ON WHICH EACH TASK IS ASIGNED
struct procName {
    int     myId;
    char*   pName;
};

/// FOR COMMUNICATING WITH MR-MPI
struct giftBox {
    int             myId;
    int             numProc;
    //int             numDbChunks;
    //int             totalDbChunks;
    char*           pName;
    //char*           queryFileName;
    //char*           dbChunkName;
    //char**          dbChunkNames;
    //int             numRound;
    //unsigned long   numQuery;
    //CBlastDemoApplication* cbDemo;
    //multimap<string, unsigned int> *mmapPName;
    //vector<string>* vPName;
    //int test;
};

/// MR-MPI CALLS
//void    get_node_name(int, KeyValue *, void *);
//void    collect_node_names(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
//void    save_node_names(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);

void    set_default_opts(CRef<CBlastOptionsHandle> optsHandle);
//void    run_blast(int itask, char *file, KeyValue *kv, void *ptr);
void    run_blast2(int itask, char *file, KeyValue *kv, void *ptr);
int     key_compare_str(char *p1, int len1, char *p2, int len2);
int     key_compare_str2(char *p1, int len1, char *p2, int len2);
int     key_compare_int(char *p1, int len1, char *p2, int len2);
bool    check_exclusion(string qGi, string sGi, int qCutLocStart, unsigned int qCutLocEnd, unsigned int sStart, unsigned int sEnd, unsigned int threshold);
void    collect_same_db(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
    
/// TOKENIZER ROUTINES
vector<string>  &split(const string &s, char delim, vector<string> &vecElems);
vector<string>  split(const string &s, char delim);

/// NOTE: NCBI PROVIDES THESE UTILS
//string        int2str(int number);
unsigned int    str2uint(string str);
int             str2int(string str);

 
int main(int argc, char **argv)
{
    /// MPI
    MPI_Init(&argc, &argv);

    char MPI_procName[MAXSTR];
    int MPI_myId, MPI_nProcs, MPI_length;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myId);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_nProcs);
    MPI_Get_processor_name(MPI_procName, &MPI_length);
    fprintf(stdout, "### [Rank %d]: %s \n", MPI_myId, MPI_procName);

    if (argc < 3) {
        if (MPI_myId == 0) printf("Syntax: mpirun -np n mrblast masterFileName outFileName [numQueryFile]\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /*
    ///MR-MPI for node name collecting
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
    //for (unsigned int i = 0; i < vecProcName.size(); ++i)
    //cout << "### Node " << MPI_myId << " name = " << vecProcName[i] << endl;
    //}

    delete mr;
    */
    
    /*
    ///
    /// Blast search v1.
    ///
    MapReduce *mr2 = new MapReduce(MPI_COMM_WORLD);
    mr2->verbosity = 0;
    mr2->timer = 1;
    mr2->mapstyle = 2;  /// master/slave mode, NOTE: the master does not run map().

    int nvecRes;
    char *masterFileName = argv[1];

    giftBox gf;
    gf.myId = MPI_myId;    
    gf.pName = MPI_procName;
    gf.dbChunkName = argv[2];
    //gf.numProc = MPI_nProcs;
    //gf.mmapPName = &mmapProcNames;
    //gf.vPName = &vecProcName;
    //gf.cbDemo = &demo;
    gf.test = 0;
    
    nvecRes = mr2->map(masterFileName, &run_blast, &gf);
    cout << "### [Rank " << MPI_myId << "] nvecRes = " << nvecRes << endl;
    cout << "gf.test = " << gf.test << endl;
    //mr2->print(-1, 1, 5, 5);
    
    //mr2->collate(NULL); //////////////////////////////////////////////
    //mr2->print(-1, 1, 5, 5);

    //nvecRes = mr2->map(mr2, &output, &gf);

    MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////

    mr2->gather(1);
    //mr2->print(0, 1, 5, 5);

    mr2->sort_keys(&key_compare_int);
    mr2->print(0, 1, 5, 5);
        
    if (MPI_myId == 0) {
        cout << "Saving...\n";        
        string outFileName = string(argv[3]) + ".txt";
        cout << "outFileName = " << outFileName << endl;
        FILE *outFile = fopen(outFileName.c_str(), "w");
        if (outFile) {
            mr2->kv->print2file(1, 5, 5, outFile);
        }
        cout << "DONE!\n";
    }

    delete mr2;
    */
    
    ///*
    ///
    /// Blast search v2
    ///
    MapReduce *mr2 = new MapReduce(MPI_COMM_WORLD);
    mr2->verbosity = 0;
    mr2->timer = 1;
    mr2->mapstyle = 2; /// master/slave mode. NOTE: the master does not run map().

    int nvecRes;
    const char *masterFileName = argv[1];

    giftBox gf;
    gf.myId = MPI_myId;    
    gf.pName = MPI_procName;
    
    ///
    /// Make a fake file for map() which contains a list of fake file 
    /// names. A file name is a form of "queryFile-dbChunkName"
    /// In map(), the file name is splitted into "queryFile" and 
    /// "dbChunkName". We've got 109 DB chunks, FYI.
    ///    
    vector<string> vWorkItems;
    if (MPI_myId == 0) {
        const char *dbChunkNameFileName = "dbchunks.txt";
        ifstream dbChunkNameFile(dbChunkNameFileName);
        string line;
        vector<string> vDbChunkNames;
        while (!getline(dbChunkNameFile, line).eof()) {
            vDbChunkNames.push_back(line);
        }
        dbChunkNameFile.close();
        unsigned int numDbChunks = vDbChunkNames.size();
        //cout << "numDbChunks = " << numDbChunks << endl;    
        ifstream masterFile(masterFileName);
        vector<string> vQueryFileNames;
        while (!getline(masterFile, line).eof()) {
            vQueryFileNames.push_back(line);
        }    
        masterFile.close();
        unsigned int numQueryFiles = vQueryFileNames.size();
        //cout << "numQueryFiles = " << numQueryFiles << endl;        
        //ofstream workItemFile(workItemFileName);            
        //if (workItemFile.is_open())     
            //for (unsigned int i = 0; i < numQueryFiles; ++i)
                //for (unsigned int j = 0; j < numDbChunks; ++j)
            ////for (unsigned int i = 0; i < 2; ++i)
                ////for (unsigned int j = 0; j < 2; ++j)
                    //workItemFile << vQueryFileNames[i] << "," << vDbChunkNames[j] << endl;
        //workItemFile.close();
        for (unsigned int i = 0; i < numQueryFiles; ++i)
            for (unsigned int j = 0; j < numDbChunks; ++j)
                vWorkItems.push_back(vQueryFileNames[i] + "," + vDbChunkNames[j]);
        cout << "### INFO: Total number of work items = " << vWorkItems.size() << endl;
        vQueryFileNames.clear();
        vDbChunkNames.clear();        
    }
    MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////
    
    char* workItemFileName = "workitems.txt";
    if (argc == 4) {        
        ///
        /// Divide the whole workload into sub jobs.
        /// Each iteration consists of calling blast and saving output
        /// into file.
        ///
        unsigned int numIters = vWorkItems.size() / str2uint(string(argv[3]));
        unsigned int numRemains = vWorkItems.size() % str2uint(string(argv[3]));
        unsigned int k = 0;
        
        for (unsigned int i = 0; i < numIters; ++i) {
            ofstream workItemFile(workItemFileName);
            if (workItemFile.is_open()) {
                for (unsigned int j = 0; j < str2uint(string(argv[3])); ++j, ++k) 
                    workItemFile << vWorkItems[k] << endl;
                workItemFile.close();
            
                //nvecRes = mr2->map(workItemFileName, &run_blast2, &gf);
                //cout << "### [Rank " << MPI_myId << "] Num res from map() = " << nvecRes << endl;
                //mr2->print(-1, 1, 5, 5);
                //mr2->collate(NULL); //////////////////////////////////////////////
                //mr2->print(-1, 1, 5, 5);
            }
        }    
        
        ofstream workItemFile2(workItemFileName);
        if (workItemFile2.is_open()) {
            for (unsigned int j = 0; j < numRemains; ++j, ++k)
                workItemFile2 << vWorkItems[k] << endl;\
            workItemFile2.close();

                //nvecRes = mr2->map(workItemFileName, &run_blast2, &gf);
                //cout << "### [Rank " << MPI_myId << "] Num res from map() = " << nvecRes << endl;
                //mr2->print(-1, 1, 5, 5);
                //mr2->collate(NULL); //////////////////////////////////////////////
                //mr2->print(-1, 1, 5, 5);            
        }        
    } 
    else { /// Process whole work items 
        if (MPI_myId == 0) {
            ofstream workItemFile(workItemFileName);
            if (workItemFile.is_open()) {
                for (unsigned int i = 0; i < vWorkItems.size(); ++i)
                    workItemFile << vWorkItems[i] << endl;
                workItemFile.close();
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////

        /// DEBUG
        char* workItemFileName2 = "test_workitems.txt";
        nvecRes = mr2->map(workItemFileName2, &run_blast2, &gf);
        ///
        //nvecRes = mr2->map(workItemFileName, &run_blast2, &gf);
        cout << "### [Rank " << MPI_myId << "] Num res from map() = " << nvecRes << endl;
        //mr2->print(-1, 1, 5, 5);
        
        mr2->sort_keys(&key_compare_str2);
        //mr2->print(-1, 1, 5, 5);
        
        //mr2->collate(NULL); //////////////////////////////////////////////
        //mr2->print(-1, 1, 5, 5);
        
        //mr2->sort_keys(&key_compare_str2);
        //mr2->print(-1, 1, 5, 5);
        
        //int nunique = mr2->reduce(&collect_same_db,NULL); /////////////////////////////////
        //mr2->print(-1, 1, 5, 5);
        
        MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////

        mr2->gather(1);
        //mr2->print(0, 1, 5, 5);
        
        mr2->sort_keys(&key_compare_str2);
        //mr2->print(0, 1, 5, 5);
        
        //nvecRes = mr2->map(mr2, &output, &gf);
        //MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////

        //mr2->gather(1);
        ////mr2->print(0, 1, 5, 5);

        //mr2->sort_keys(&key_compare_int);
        //mr2->print(0, 1, 5, 5);
            
        if (MPI_myId == 0) {
            cout << "Saving...\n";        
            string outFileName = string(argv[2]) + ".txt";
            cout << "outFileName = " << outFileName << endl;
            FILE *outFile = fopen(outFileName.c_str(), "w");
            if (outFile) {
                mr2->kv->print2file(1, 5, 5, outFile);
            }
            cout << "DONE!\n";
        }
    }

    delete mr2;
    //*/
    
    MPI::Finalize();
}

/* ---------------------------------------------------------------------- */
void run_blast2(int itask, char *file, KeyValue *kv, void *ptr)
/* ---------------------------------------------------------------------- */
{
    giftBox *gf = (giftBox *) ptr;
    
    ///
    /// MAKE A OPTION HANDLE AND CBLASTINPUTSOURCE
    ///
    EProgram program = ProgramNameToEnum("blastn");
    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));
    set_default_opts(opts);
    opts->Validate();

    /// DEBUG
    //if (gf->myId == 0) 
        //opts->GetOptions().DebugDumpText(cerr, "opts", 1);

    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
        throw std::runtime_error("Could not initialize object manager");
    }

    bool isProtein = false;
    SDataLoaderConfig dlconfig(isProtein);
    CBlastInputSourceConfig iconfig(dlconfig);
    fprintf(stdout, "### INFO: [Rank %d]: %s, itask = %d, Work item = %s\n",
            gf->myId, gf->pName, itask, file);
    //string qFileName = "test2.query";
    //ifstream queryFile(qFileName.c_str());
    //CBlastFastaInputSource fasta_input(qFileName, iconfig);
    
    /// Split 
    string workItem(file);
    vector<string> vWorkItem  = split(workItem, ',');
    ifstream queryFile(vWorkItem[0].c_str());

    /// USE QUERY STRING INSTEAD OF FILE INPUT
    string header, seq;
    vector<string> vecSeq;

    ///
    /// TARGET DB CHUNK SETTING
    ///
    string dbChunkName = vWorkItem[1];
    const CSearchDatabase target_db(dbChunkName, CSearchDatabase::eBlastDbIsNucleotide);
    
    ///
    /// READ SEQ(S) AND RUN BLAST
    ///
    while (!getline(queryFile, header).eof()) {
        getline(queryFile, seq);
        
        vector<string> vHeaders;
        
        if (header.length() == 0 || seq.length() == 0) {
            fprintf(stderr, "### ERROR: [Rank %d]: %s, itask = %d, \
                file = %s, Seq read error.\n",
                gf->myId, gf->pName, itask, file);
            exit(0);
        }
        //fprintf(stdout, "### INFO: [Rank %d]: %s, itask = %d, %s\n",
            //gf->myId, gf->pName, itask, header.c_str());

        string query = header + '\n' + seq;
        vHeaders.push_back(header);
        
        ///
        /// ACCUMMULATE SEQS WHEN MAXQUERYNUM > 1
        ///
        unsigned int j = 1;
        while (j < MAXQUERYNUM && !getline(queryFile, header).eof()) {
            getline(queryFile, seq);
            //fprintf(stdout, "### INFO: [Rank %d]: %s, itask = %d, %s, num queries = %d\n",
                //gf->myId, gf->pName, itask, header.c_str(), j + 1);
            query = query + '\n' + header + '\n' + seq;
            vHeaders.push_back(header);
            j++;
        }
        
        /// TOKENIZE QUERY HEADER
        //vector<string> vQeuryId  = split(header, '|');
        //string gi                = vQeuryId[1];             /// GI
        //string qid               = vQeuryId[2];             /// QUERY ID
        //unsigned int oLen              = str2uint(vQeuryId[3]);   /// LENGTH OF THE ORIG SEQ
        //int qStart               = str2int(vQeuryId[4]);    /// CUT COORDINATES - START
        //unsigned int qEnd              = str2uint(vQeuryId[5]);   /// CUT COORDINATES - END
        //fprintf(stderr, "### INFO: [Rank %d]: %s, itask = %d, %s %s %d %d %d\n",
                //gf->myId, gf->pName, itask, gi.c_str(), qid.c_str(), oLen, qStart, qEnd);
                        
        ///
        /// SET QUERIES AS FASTA INPUT
        ///
        CBlastFastaInputSource fasta_input(query, iconfig);
        CBlastInput blastInput(&fasta_input);
        CScope scope(*objmgr);
        TSeqLocVector queryLoc = blastInput.GetAllSeqLocs(scope);
        CRef<IQueryFactory> queryFactory(new CObjMgr_QueryFactory(queryLoc));

        ///
        /// RUN BLAST
        ///
        CLocalBlast blaster(queryFactory, opts, target_db);
        CSearchResultSet results = *blaster.Run();

        ///
        /// GET THE RESULTS
        ///

        /// GET WARNING MESSAGES.
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
        //cout << "results.GetNumQueries() = " << results.GetNumQueries() << endl;
        for (unsigned int i = 0; i < results.GetNumResults(); i++) {

            //CConstRef<CSeq_id> seq_id = results[i].GetSeqId();
            //cout << "qID = " << seq_id->AsFastaString() << endl;

            CConstRef<CSeq_align_set> aln_set = results[i].GetSeqAlign();
            //cout << MSerial_AsnText << *aln_set;
            //cout << MSerial_Xml << *aln_set;
            
            ///
            /// PRINT TABULAR FORMAT
            ///                                            
            if (results[i].HasAlignments()) {
                ITERATE(CSeq_align_set::Tdata, itr, aln_set->Get()) {
                    const CSeq_align& s = **itr;
                    //tabinfo.SetFields(s, scope, &m_ScoringMatrix);
                    //tabinfo.Print();
                    
                    string qNum = s.GetSeq_id(QUERY).GetSeqIdString();    
                    //cout << "qNum = " << qNum << endl;       
                    //cout << queryFactory.GetQueryInfo() << endl;
                    //cout << s.GetSeq_id(QUERY).AsFastaString() << endl;
                    string seqId = s.GetSeq_id(SUBJECT).GetSeqIdString();
                    
                    //double pIdentityGapped, pIdentityUngapped, pIdentityGapOpeningOnly, pCoverage;
                    //s.GetNamedScore(CSeq_align::eScore_PercentIdentity_Gapped, pIdentityGapped);
                    //s.GetNamedScore(CSeq_align::eScore_PercentIdentity_Ungapped, pIdentityUngapped);
                    //s.GetNamedScore(CSeq_align::eScore_PercentIdentity_GapOpeningOnly, pIdentityGapOpeningOnly);
                    //s.GetNamedScore(CSeq_align::eScore_PercentCoverage, pCoverage);
                    
                    unsigned int alignLen, qStart, qEnd, sStart, sEnd;
                    //CRange<TSeqPos> qRange, sRange;
                    //qRange = s.GetSeqRange(QUERY);
                    //sRange = s.GetSeqRange(SUBJECT);
                    qStart      = s.GetSeqStart(QUERY);
                    qEnd        = s.GetSeqStop(QUERY);
                    sStart      = s.GetSeqStart(SUBJECT);
                    sEnd        = s.GetSeqStop(SUBJECT);
                    alignLen    = s.GetAlignLength();
                    
                    double eValue;
                    int genericScore, bitScore;
                    s.GetNamedScore(CSeq_align::eScore_EValue, eValue);
                    s.GetNamedScore(CSeq_align::eScore_Score, genericScore);                    
                    s.GetNamedScore(CSeq_align::eScore_BitScore, bitScore);       
                    
                    ///
                    /// TOKENIZE QUERY HEADER
                    ///
                    string qHeader = vHeaders[str2uint(qNum)-1];
                    //cout << "qHeader = " << qHeader << endl;
                    vector<string> vQeuryId  = split(qHeader, '|');
                    string qGi               = vQeuryId[1];             /// GI
                    string qid               = vQeuryId[2];             /// QUERY ID
                    //unsigned int origLen     = str2uint(vQeuryId[3]);   /// LENGTH OF THE ORIG SEQ
                    int qCutLocStart         = str2int(vQeuryId[4]);    /// CUT COORDINATES - START
                    unsigned int qCutLocEnd  = str2uint(vQeuryId[5]);   /// CUT COORDINATES - END
                    //fprintf(stderr, "### INFO: [Rank %d]: %s, itask = %d, %s %s %d %d %d\n",
                        //gf->myId, gf->pName, itask, qGi.c_str(), qid.c_str(), 
                        //origLen, qCutLocStart, qCutLocEnd);
                        
                    ///
                    /// ADD A CSV BLAST RESULT TO KV
                    ///
                    //int cutoffScore = 100;
                    //double eValueThreshold = 10.0;
                    
                    //if (bitScore >= cutoffScore && eValue < eValueThreshold) {
                        //fprintf(stdout, "%s,%s,%s, %f,%f,%f,%f, %d,%d,%d,%f,%d\n", 
                        //fprintf(stdout, "%s,%s,%s,%d,%d,%d,%1.e,%d,%d\n", 
                            //qid.c_str(), qHeader.substr(1,1000).c_str(), seqId.c_str(), 
                            ////pIdentityGapped, pIdentityUngapped, pIdentityGapOpeningOnly, pCoverage,
                            //alignLen, qStart, qEnd, 
                            //eValue, genericScore, bitScore);
                    //}
                    if (!check_exclusion(qGi, seqId, qCutLocStart, qCutLocEnd, sStart, sEnd, 100)) {
                        
                        char blastRes[MAXSTR];
                        
                        /// DEBUG
                        //fprintf(stdout, "%s,%s,%s,%d,%d,%d,%1.e,%d\n", 
                            //qid.c_str(), qHeader.substr(1,MAXSTR).c_str(), seqId.c_str(), 
                            //alignLen, qStart, qEnd, 
                            ////eValue, genericScore, bitScore);
                            //eValue, bitScore);
                        
                        sprintf(blastRes, "%s,%s,%d,%d,%d,%1.e,%d", 
                            qHeader.substr(1,MAXSTR).c_str(), seqId.c_str(), 
                            alignLen, qStart, qEnd, 
                            //eValue, genericScore, bitScore);
                            eValue, bitScore);
                        
                        char *newKey = (char*)((qid+','+dbChunkName).c_str());

                        /// ADD <KEY = "QUERYID,CHUNKNAME", VALUE="BLASTvecResULT"> TO KV
                        kv->add(newKey, strlen(newKey) + 1, blastRes, strlen(blastRes) + 1);
                    }
                    else 
                        fprintf(stderr, "### INFO: EXCLUSION - qGi=%s sGi=%s \
                            qCutLocStart=%d, qCutLocEnd=%d, sStart=%d, sEnd=%d\n",
                            qGi.c_str(), seqId.c_str(), qCutLocStart, qCutLocEnd, sStart, sEnd);
                        
                }
            }
        }
    } /// END WHILE
}


/* ---------------------------------------------------------------------- */
//void run_blast(int itask, char *file, KeyValue *kv, void *ptr)
/* ---------------------------------------------------------------------- */
/*
{
    giftBox *gf = (giftBox *) ptr;
    
    ///
    /// MAKE A OPTION HANDLE AND CBLASTINPUTSOURCE
    ///
    EProgram program = ProgramNameToEnum("blastn");
    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));
    set_default_opts(opts);
    opts->Validate();

    /// DEBUG
    //if (gf->myId == 0) 
        //opts->GetOptions().DebugDumpText(cerr, "opts", 1);

    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
        throw std::runtime_error("Could not initialize object manager");
    }

    bool isProtein = false;
    SDataLoaderConfig dlconfig(isProtein);
    CBlastInputSourceConfig iconfig(dlconfig);
    ifstream queryFile(file);
    cout << "qFile Name = " << file << endl;
    //string qFileName = "test2.query";
    //ifstream queryFile(qFileName.c_str());
    //CBlastFastaInputSource fasta_input(qFileName, iconfig);

    /// USE QUERY STRING INSTEAD OF FILE INPUT
    string header, seq;
    vector<string> vecSeq;

    ///
    /// TARGET DB SETTING
    ///
    cout << "dbChunkName = " << gf->dbChunkName << endl;
    string dbChunkName(gf->dbChunkName);
    const CSearchDatabase target_db(dbChunkName, CSearchDatabase::eBlastDbIsNucleotide);
    
    ///
    /// READ SEQ(S) AND RUN BLAST
    ///
    while (!getline(queryFile, header).eof()) {
        getline(queryFile, seq);
        
        vector<string> vHeaders;
        
        if (header.length() == 0 || seq.length() == 0) {
            fprintf(stderr, "### ERROR: [Rank %d]: %s, itask = %d, \
                file = %s, Seq read error.\n",
                gf->myId, gf->pName, itask, file);
            exit(0);
        }
        fprintf(stdout, "### INFO: [Rank %d]: %s, itask = %d, %s\n",
            gf->myId, gf->pName, itask, header.c_str());

        string query = header + '\n' + seq;
        vHeaders.push_back(header);
        
        ///
        /// ACCUMMULATE SEQS WHEN MAXQUERYNUM > 1
        ///
        unsigned int j = 1;
        while (j < MAXQUERYNUM && !getline(queryFile, header).eof()) {
            getline(queryFile, seq);
            fprintf(stdout, "### INFO: [Rank %d]: %s, itask = %d, %s, num queries = %d\n",
                gf->myId, gf->pName, itask, header.c_str(), j + 1);
            query = query + '\n' + header + '\n' + seq;
            vHeaders.push_back(header);
            j++;
        }
        
        /// TOKENIZE QUERY HEADER
        //vector<string> vQeuryId  = split(header, '|');
        //string gi                = vQeuryId[1];             /// GI
        //string qid               = vQeuryId[2];             /// QUERY ID
        //unsigned int oLen              = str2uint(vQeuryId[3]);   /// LENGTH OF THE ORIG SEQ
        //int qStart               = str2int(vQeuryId[4]);    /// CUT COORDINATES - START
        //unsigned int qEnd              = str2uint(vQeuryId[5]);   /// CUT COORDINATES - END
        //fprintf(stderr, "### INFO: [Rank %d]: %s, itask = %d, %s %s %d %d %d\n",
                //gf->myId, gf->pName, itask, gi.c_str(), qid.c_str(), oLen, qStart, qEnd);
                        
        ///
        /// SET QUERIES AS FASTA INPUT
        ///
        CBlastFastaInputSource fasta_input(query, iconfig);
        CBlastInput blastInput(&fasta_input);
        CScope scope(*objmgr);
        TSeqLocVector queryLoc = blastInput.GetAllSeqLocs(scope);
        CRef<IQueryFactory> queryFactory(new CObjMgr_QueryFactory(queryLoc));

        ///
        /// RUN BLAST
        ///
        CLocalBlast blaster(queryFactory, opts, target_db);
        CSearchResultSet results = *blaster.Run();

        ///
        /// GET THE RESULTS
        ///

        /// GET WARNING MESSAGES.
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
        cout << "results.GetNumQueries() = " << results.GetNumQueries() << endl;
        for (unsigned int i = 0; i < results.GetNumResults(); i++) {

            //CConstRef<CSeq_id> seq_id = results[i].GetSeqId();
            //cout << "qID = " << seq_id->AsFastaString() << endl;

            CConstRef<CSeq_align_set> aln_set = results[i].GetSeqAlign();
            //cout << MSerial_AsnText << *aln_set;
            //cout << MSerial_Xml << *aln_set;
            
            ///
            /// PRINT TABULAR FORMAT
            ///                                            
            if (results[i].HasAlignments()) {
                ITERATE(CSeq_align_set::Tdata, itr, aln_set->Get()) {
                    const CSeq_align& s = **itr;
                    //tabinfo.SetFields(s, scope, &m_ScoringMatrix);
                    //tabinfo.Print();
                    
                    string qNum = s.GetSeq_id(QUERY).GetSeqIdString();    
                    //cout << "qNum = " << qNum << endl;       
                    //cout << queryFactory.GetQueryInfo() << endl;
                    //cout << s.GetSeq_id(QUERY).AsFastaString() << endl;
                    string seqId = s.GetSeq_id(SUBJECT).GetSeqIdString();
                    
                    //double pIdentityGapped, pIdentityUngapped, pIdentityGapOpeningOnly, pCoverage;
                    //s.GetNamedScore(CSeq_align::eScore_PercentIdentity_Gapped, pIdentityGapped);
                    //s.GetNamedScore(CSeq_align::eScore_PercentIdentity_Ungapped, pIdentityUngapped);
                    //s.GetNamedScore(CSeq_align::eScore_PercentIdentity_GapOpeningOnly, pIdentityGapOpeningOnly);
                    //s.GetNamedScore(CSeq_align::eScore_PercentCoverage, pCoverage);
                    
                    unsigned int alignLen, qStart, qEnd, sStart, sEnd;
                    //CRange<TSeqPos> qRange, sRange;
                    //qRange = s.GetSeqRange(QUERY);
                    //sRange = s.GetSeqRange(SUBJECT);
                    qStart      = s.GetSeqStart(QUERY);
                    qEnd        = s.GetSeqStop(QUERY);
                    sStart      = s.GetSeqStart(SUBJECT);
                    sEnd        = s.GetSeqStop(SUBJECT);
                    alignLen    = s.GetAlignLength();
                    
                    double eValue;
                    int genericScore, bitScore;
                    s.GetNamedScore(CSeq_align::eScore_EValue, eValue);
                    s.GetNamedScore(CSeq_align::eScore_Score, genericScore);                    
                    s.GetNamedScore(CSeq_align::eScore_BitScore, bitScore);       
                    
                    ///
                    /// TOKENIZE QUERY HEADER
                    ///
                    string qHeader = vHeaders[str2uint(qNum)-1];
                    //cout << "qHeader = " << qHeader << endl;
                    vector<string> vQeuryId  = split(qHeader, '|');
                    string qGi               = vQeuryId[1];             /// GI
                    string qid               = vQeuryId[2];             /// QUERY ID
                    unsigned int origLen           = str2uint(vQeuryId[3]);   /// LENGTH OF THE ORIG SEQ
                    int qCutLocStart         = str2int(vQeuryId[4]);    /// CUT COORDINATES - START
                    unsigned int qCutLocEnd        = str2uint(vQeuryId[5]);   /// CUT COORDINATES - END
                    fprintf(stderr, "### INFO: [Rank %d]: %s, itask = %d, %s %s %d %d %d\n",
                        gf->myId, gf->pName, itask, qGi.c_str(), qid.c_str(), 
                        origLen, qCutLocStart, qCutLocEnd);
                        
                    ///
                    /// ADD A CSV BLAST RESULT TO KV
                    ///
                    //int cutoffScore = 100;
                    //double eValueThreshold = 10.0;
                    
                    //if (bitScore >= cutoffScore && eValue < eValueThreshold) {
                        //fprintf(stdout, "%s,%s,%s, %f,%f,%f,%f, %d,%d,%d,%f,%d\n", 
                        //fprintf(stdout, "%s,%s,%s,%d,%d,%d,%1.e,%d,%d\n", 
                            //qid.c_str(), qHeader.substr(1,1000).c_str(), seqId.c_str(), 
                            ////pIdentityGapped, pIdentityUngapped, pIdentityGapOpeningOnly, pCoverage,
                            //alignLen, qStart, qEnd, 
                            //eValue, genericScore, bitScore);
                    //}
                    if (!check_exclusion(qGi, seqId, qCutLocStart, qCutLocEnd, sStart, sEnd, 100)) {
                        
                        char blastRes[MAXSTR];
                        
                        /// DEBUG
                        fprintf(stdout, "%s,%s,%s,%d,%d,%d,%1.e,%d,%d\n", 
                            qid.c_str(), qHeader.substr(1,MAXSTR).c_str(), seqId.c_str(), 
                            alignLen, qStart, qEnd, 
                            eValue, genericScore, bitScore);
                        
                        sprintf(blastRes, "%s,%s,%d,%d,%d,%1.e,%d,%d", 
                            qHeader.substr(1,MAXSTR).c_str(), seqId.c_str(), 
                            alignLen, qStart, qEnd, 
                            eValue, genericScore, bitScore);
                        
                        char *newKey = (char*)(qid.c_str());

                        /// ADD <KEY = "QUERYID,CHUNKNAME", VALUE="BLASTvecResULT"> TO KV
                        kv->add(newKey, strlen(newKey) + 1, blastRes, strlen(blastRes) + 1);
                    }
                    else 
                        fprintf(stderr, "### INFO: EXCLUSION - qGi=%s sGi=%s \
                            qCutLocStart=%d, qCutLocEnd=%d, sStart=%d, sEnd=%d\n",
                            qGi.c_str(), seqId.c_str(), qCutLocStart, qCutLocEnd, sStart, sEnd);
                        
                }
            }
        }
    } // END WHILE
}
*/

/* ---------------------------------------------------------------------- */
void collect_same_db(char *key, int keybytes, char *multivalue,
    int nvalues, int *valuebytes, KeyValue *kv, void *ptr) 
/* ---------------------------------------------------------------------- */
{
    //kv->add(key,keybytes,(char *) &nvalues,sizeof(int));
    unsigned int sum = 0;
    for (int i = 0; i < nvalues; ++i)
        sum += valuebytes[i];    
    ///
    /// Why?
    /// multivalue has null characters inside
    ///
    //char temp[sum];
    //cout << sum << endl;
    //memcpy(temp, multivalue, sum);
    //int idx = 0;
    //for (int i = 0; i < nvalues; ++i) {
        //char temp2[sum];
        //memcpy(temp, multivalue, valuebytes[i]);
        //idx += valuebytes[i];
        //cout << temp2 << endl;
    //}
    //cout << endl;
    //kv->add(key,keybytes,temp,sum);
    kv->add(key,keybytes,multivalue,sum);
}

/* ---------------------------------------------------------------------- */
bool check_exclusion(string qGi, string sGi, int qCutLocStart, 
    unsigned int qCutLocEnd, unsigned int sStart, unsigned int sEnd, 
    unsigned int threshold)
/* ---------------------------------------------------------------------- */
{   
    ///
    /// To exclude Blast result from the original sequence from which 
    /// the input query is originated (sampled). Basically if qGi == sGi 
    /// and qCutLocStart is similar with sStart and qCutLocEnd is similar 
    /// with sEnd in terms of coordinates, the result should be excluded.
    /// 
    /// Orig seq: ----------------XXXXXXXXXXXXXXXX----------------------
    ///                           |              |
    ///                    qCutLocStart      qCutLocEnd 
    ///
    /// Query:                    XXXXXXXXXXXXXXXX
    ///                             |          | 
    ///                          qStart       qEnd
    ///
    /// Subject:  ------------------XXXXXXXXXXXX------------------------
    ///                             |          |
    ///                           sStart      sEnd
    ///
    bool ret = false;
    
    if (qGi == sGi) {
        if (qCutLocStart < 0) {            
            /// In >gi|222299657|18|3605|-400|3604
            /// -400|3604 means query[-400:3604] in Python.
            qCutLocStart = qCutLocEnd + 1 - qCutLocStart;
            qCutLocEnd += 1;
        }        
        if ((qCutLocStart - threshold <= sStart && sStart <= qCutLocStart + threshold) &&
            (qCutLocEnd - threshold <= sEnd && sEnd <= qCutLocEnd + threshold)) 
            ret = true;
    }    
    
    return ret;
}

/* ---------------------------------------------------------------------- */
void set_default_opts(CRef<CBlastOptionsHandle> optsHandle)
/* ---------------------------------------------------------------------- */
{
    optsHandle->SetEvalueThreshold(10);
    //optsHandle->SetMatchReward(0);
    //optsHandle->SetMismatchPenalty(0);
    //optsHandle->SetMatrixName("BLOSUM62");
    optsHandle->SetCutoffScore(50);
    if (CBlastNucleotideOptionsHandle* nuclHandle =
            dynamic_cast<CBlastNucleotideOptionsHandle*>(&*optsHandle)) {

        //nuclHandle->SetMatchReward(0);
        //nuclHandle->SetMismatchPenalty(0);
        nuclHandle->SetMatrixName("BLOSUM62");
        nuclHandle->SetCutoffScore(50);
        nuclHandle->SetEvalueThreshold(10);
    }

    return;
}

/* ---------------------------------------------------------------------- */
int key_compare_str(char *p1, int len1, char *p2, int len2)
/* ---------------------------------------------------------------------- */
{
    return strcmp(p1, p2);
}

/* ---------------------------------------------------------------------- */
int key_compare_str2(char *p1, int len1, char *p2, int len2)
/* ---------------------------------------------------------------------- */
{
    string str1(p1);
    string str2(p2);
    vector<string> vStr1 = split(str1, ',');
    vector<string> vStr2 = split(str2, ',');
    unsigned int i1 = str2uint(vStr1[0]);
    unsigned int i2 = str2uint(vStr2[0]);
    //cout << "i1 vs i2 = " << i1 << " " << i2 << endl;
    //return strcmp(p1, p2);
    /// Ascending order
    if (i1 > i2) return 1;
    else if (i1 < i2) return -1;
    /// Descending order
    //if (i1 > i2) return -1;
    //else if (i1 < i2) return 1;
    
    //else return 0;
    else vStr1[1].compare(vStr2[1]);
    
}

/* ---------------------------------------------------------------------- */
int key_compare_int(char *p1, int len1, char *p2, int len2)
/* ---------------------------------------------------------------------- */
{
    //int i1 = *(int *) p1;
    //int i2 = *(int *) p2;
    unsigned int i1 = str2uint(string(p1));
    unsigned int i2 = str2uint(string(p2));
    /// Ascending order
    if (i1 > i2) return 1;
    else if (i1 < i2) return -1;
    /// Descending order
    //if (i1 > i2) return -1;
    //else if (i1 < i2) return 1;
    else return 0;
}

///* ---------------------------------------------------------------------- */
//void get_node_name(int itask, KeyValue *kv, void *ptr)
///* ---------------------------------------------------------------------- */
//{
    //procName *nn = (procName *) ptr;
    //char id[MAXSTR];
    //sprintf(id, "%d", nn->myId);
    //kv->add(nn->pName, strlen(nn->pName) + 1, id, strlen(id) + 1);
//}

/* ---------------------------------------------------------------------- */
//void collect_node_names(char *key, int keybytes, char *multivalue,
//                        int nvalues, int *valuebytes, KeyValue *kv,
//                        void *ptr)
/* ---------------------------------------------------------------------- */
//{
    //kv->add(key, keybytes, multivalue, sizeof(multivalue));
//}

///* ---------------------------------------------------------------------- */
//void save_node_names(uint64_t itask, char *key, int keybytes, char *value,
                     //int valuebytes, KeyValue *kv, void *ptr)
///* ---------------------------------------------------------------------- */
//{
    ////multimap<string, unsigned int> *nNames = (multimap<string, unsigned int> *) ptr;
    //vector<string>* vecNames = (vector<string> *) ptr;
    //string pName(key);
    ////unsigned int myid = atoi(value);
    ////nNames->insert(pair<string, unsigned int>(pName, myid));
    //vector<string>::const_iterator loc = find(vecNames->begin(), vecNames->end(), pName);
    //if (loc == vecNames->end())
        //vecNames->push_back(pName);
//}

/// TOKENIZER ROUTINES
vector<string> &split(const string &s, char delim, vector<string> &vecElems)
{
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        vecElems.push_back(item);
    }
    return vecElems;
}

vector<string> split(const string &s, char delim)
{
    vector<string> vecElems;
    return split(s, delim, vecElems);
}

/// NOTE: NCBI PROVIDES THESE UTILS
/// INT TO STRING 
//string int2str(int number)
//{
    //stringstream ss;
    //ss << number;
    //return ss.str();
//}

unsigned int str2uint(string str)
{
    std::stringstream ss;
    ss << str;
    unsigned int f;
    ss >> f;
    return f;
}

int str2int(string str)
{
    std::stringstream ss;
    ss << str;
    int f;
    ss >> f;
    return f;
}

/// EOF

