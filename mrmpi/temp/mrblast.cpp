////////////////////////////////////////////////////////////////////////
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
//
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
//      10.01.2010  Value sorting by bit score is done.
//
//      10.06.2010  Divide & conquer for KV overflow is done.
//
//      10.11.2010  Instrumenting for section timing.
//
//      10.13.2010  Doxygen commenting
//                  Add map() in mapreduce.cpp for handling vector<string> 
//                  work items.
//
//      10.20.2010  Save in workers.
//
//      10.27.2010  Custom scheduler done!
//
//
////////////////////////////////////////////////////////////////////////

/*
 *     This program is free software; you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation; either version 2 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program; if not, write to the Free Software
 *     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *     MA 02110-1301, USA.
 */

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

/// ASN & split
#include <sstream>

/// blastn importing
//#include <algo/blast/blastinput/blastn_args.hpp>
//#include <algo/blast/format/blast_format.hpp>
//#include "blast_app_util.hpp"

/// Tabular
//#include <algo/blast/format/blast_format.hpp>
//#include <algo/blast/format/blastxml_format.hpp>
//#include <algo/blast/format/data4xmlformat.hpp> /// NCBI_FAKE_WARNING
//#include <algo/blast/format/build_archive.hpp>
//#include <objects/seq/Seq_annot.hpp>
//#include <objects/general/User_object.hpp>
//#include <objects/general/User_field.hpp>
//#include <algo/blast/nodeIdx/blast_stat.h>
//#include <corelib/ncbiutil.hpp>                 /// for FindBestChoice

/// XML
//#include <misc/xmlwrapp/xmlwrapp.hpp>           /// for XML::
//#include <misc/xmlwrapp/xsltwrapp.hpp>

///For CSeq_align for output
#include <objects/seqalign/Seq_align.hpp>
#include <util/range.hpp>

/// For typedef unsigned long long int uint64_t
#include <stdint.h>

/// For timing
#include <sys/time.h>
#define DEBUG 1

/// For queue
#include <queue>

USING_NCBI_SCOPE;
USING_SCOPE(blast);
//USING_SCOPE(objects);
//USING_SCOPE(align_format);
//USING_SCOPE(xml);

const uint64_t MAXSTR = 256;
const uint64_t MAXQUERYNUM = 10000;  /// num query to accumulate from qeury file
const uint64_t QUERY = 0;
const uint64_t SUBJECT = 1;
const uint64_t EXCLUSION = 100; /// Exclusion threshold = 100bp

const uint64_t CUTOFFSCORE = 100;
const float EVALUE = 1.0;

const int NUMTOTALDBCHUNKS = 109;

/// TO COLLECT PROC NAMES ON WHICH EACH TASK IS ASIGNED
struct procName {
    int myId;
    char* pName;
};

/// FOR COMMUNICATING WITH MR-MPI
struct giftBox {
    int myId;
    char* pName;
    //int numProc;    
    //char* outFileName;
    //string workItemProcessed;
    //vector<string> *dbName;
    //multimap<string, int>  *mDbCoreNum;
    //vector<string> *vDbChunkNames;
    //queue<string> *qWorkItems;
};

/// MR-MPI CALLS
//void    get_node_name(int, KeyValue *, void *);
//void    collect_node_names(char *key, int keybytes, char *multivalue, 
                             //int nvalues, int *valuebytes, KeyValue *kv, 
                             //void *ptr);
//void    save_node_names(uint64_t itask, char *key, int keybytes, char *value, 
                          //int valuebytes, KeyValue *kv, void *ptr);

void    set_default_opts(CRef<CBlastOptionsHandle> optsHandle);
//void    run_blast(int itask, char *file, KeyValue *kv, void *ptr);
void    run_blast2(int itask, char *file, KeyValue *kv, void *ptr);
//int     key_compare_str(char *p1, int len1, char *p2, int len2);
//int     key_compare_str2(char *p1, int len1, char *p2, int len2);
//int     key_compare_int(char *p1, int len1, char *p2, int len2);
bool    check_exclusion(string qGi, string sGi, int qCutLocStart, 
                        uint64_t qCutLocEnd, uint64_t sStart, uint64_t sEnd, 
                        uint64_t threshold);

void    remove_DBname_from_key(char *key, int keybytes, char *multivalue, int nvalues, 
                        int *valuebytes, KeyValue *kv, void *ptr);
//void    sort_multivalues_by_score(char *key, int keybytes, char *multivalue, int nvalues, 
                        //int *valuebytes, KeyValue *kv, void *ptr);
void    sort_multivalues_by_score2(char *key, int keybytes, char *multivalue, int nvalues, 
                        int *valuebytes, KeyValue *kv, void *ptr);                      
//void    save_output(uint64_t itask, char *key, int keybytes, char *value,
                        //int valuebytes, KeyValue *kv, void *ptr);

//void    AtExit(void);
        
/// TOKENIZER ROUTINES
vector<string>  &split(const string &s, char delim, vector<string> &vecElems);
vector<string>  split(const string &s, char delim);

/// NOTE: NCBI PROVIDES THESE UTILS
//string        int2str(int number);
string          uint2str(uint64_t number);
uint64_t        str2uint(string str);
int             str2int(string str);

/* ---------------------------------------------------------------------- */
int main(int argc, char **argv)
/* ---------------------------------------------------------------------- */
{
    /// To ensure MPI::Finalize();
    //atexit(AtExit);

    /// Declare time variables:
    #ifdef DEBUG
    time_t  t0, t1, t2, t3, accMapWTime=0; 
    t0 = time(NULL);
    clock_t c0, c1, c2, c3, accMapCTime=0;     
    c0 = clock();
    #endif

    /// MPI
    MPI_Init(&argc, &argv);

    char MPI_procName[MAXSTR];
    int MPI_myId, MPI_nProcs, MPI_length;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myId);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_nProcs);
    MPI_Get_processor_name(MPI_procName, &MPI_length);
    fprintf(stdout, "### INFO: [Rank %d] %s \n", MPI_myId, MPI_procName);
    MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////
   
//try{
    ///
    /// If numQueryFile = 0, whole work item file are processed.
    /// If numQeuryFile > 1, work item file are splitted into sub work item lists
    /// and the sub work item files are processed one by one.
    /// If mapstyle = 2, to use the default master/slave mode in MR-MPI.
    /// If mapstyle = 3, to use the modified master/slave mode for controlling 
    /// DB chunk locality.
    ///
    if (argc < 3) {
        if (MPI_myId == 0) printf("Syntax: mpirun -np n mrblast \
            masterFileName FileNamePrefix numQueryFile=0 mapstyle=2\n");
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

    //multimap<string, uint64_t> mmapProcNames;
    vector<string> vecProcName;
    mr->map(mr, &save_node_names, &vecProcName);
    //cout << vecProcName.size() << endl;

    MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////

    /// DEBUG
    //if (MPI_myId == 0) {
    ////pair<multimap<string, uint64_t>::iterator, multimap<string, uint64_t>::iterator> ppp;
    ////ppp = mmapProcNames.equal_range("ssul-lx");
    ////for (multimap<string, uint64_t>::iterator it2 = ppp.first; it2 != ppp.second; ++it2) {
    ////cout << "  [" << (*it2).first << ", " << (*it2).second << "]" << endl;
    ////}
    //for (uint64_t i = 0; i < vecProcName.size(); ++i)
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

    /// 
    /// Blast search v2
    ///
    MapReduce *mr2 = new MapReduce(MPI_COMM_WORLD);
    mr2->verbosity = 0;
    mr2->timer = 0;
    mr2->mapstyle = atoi(argv[4]);  /// master/slave mode.
    MPI_Barrier(MPI_COMM_WORLD);

    uint64_t nvecRes;
    const char *masterFileName = argv[1];
    giftBox gf;
    gf.myId = MPI_myId;
    gf.pName = MPI_procName;
    string prefix(argv[2]);

    ///
    /// Make a fake file for map() which contains a list of fake file
    /// names. A file name is a form of "queryFile,dbChunkName"
    /// In map(), the file name is splitted into "queryFile" and
    /// "dbChunkName". We've got 109 DB chunks, FYI.
    ///
    vector<string> vWorkItems;
    vector<string> vDbChunkNames;
    vector<string> vQueryFileNames;
    uint64_t numDbChunks;
    uint64_t numQueryFiles;
    //queue<string> qWorkItems;
    
    if (MPI_myId == 0) {
        /// DEBUG
        //const char *dbChunkNameFileName = "test_dbchunks.txt";
        ///
        const char *dbChunkNameFileName = "dbchunks.txt";
        ifstream dbChunkNameFile(dbChunkNameFileName);
        if (!dbChunkNameFile.is_open()) {
            cerr << "### ERROR: dbchunks.txt open error.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        string line;
        //vector<string> vDbChunkNames;
        while (!getline(dbChunkNameFile, line).eof()) {
            vDbChunkNames.push_back(line);
        }
        dbChunkNameFile.close();        
        numDbChunks = vDbChunkNames.size();
        
        ifstream masterFile(masterFileName);
        if (!masterFile.is_open()) {
            cerr << "### ERROR: masterFile open error.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        //vector<string> vQueryFileNames;
        while (!getline(masterFile, line).eof()) {
            vQueryFileNames.push_back(line);
        }
        masterFile.close();
        numQueryFiles = vQueryFileNames.size();
        for (uint64_t i = 0; i < numQueryFiles; ++i) {
            for (uint64_t j = 0; j < numDbChunks; ++j) { 
                vWorkItems.push_back(vQueryFileNames[i] + "," + vDbChunkNames[j]);
                //qWorkItems.push(vQueryFileNames[i] + "," + vDbChunkNames[j]);
            }
        }
        cout << "### INFO: Total number of work items = "
             << vWorkItems.size() << endl;
        vQueryFileNames.clear();
        vDbChunkNames.clear();
        //gf.vDbChunkNames = &vDbChunkNames;
        //gf.qWorkItems = &qWorkItems;
        
        /*
        ///
        /// Prepare a mutimap for mapstyle=3 in mapreduce.cpp
        ///
        int numNodes = MPI_nProcs / 4;
        int numDBPerNode = int(ceil(NUMTOTALDBCHUNKS / numNodes));
        int numNodeRemain = MPI_nProcs - NUMTOTALDBCHUNKS;
        multimap<string, int> mDbCoreNum;
        int nodeIdx = 0;
        int cnt = 0;
        for (uint64_t i = 0; i < vDbChunkNames.size(); ++i) {
            //mDbCoreNum[vDbChunkNames[i]] = nodeIdx;
            mDbCoreNum.insert(pair <string, int> (vDbChunkNames[i], nodeIdx)); 
            if (cnt == numDBPerNode - 1) {
                nodeIdx++;
                cnt = 0;
            }
            else cnt++;
        }
        
        /// 
        /// if nproc > number of DB chunks, just assign from the start.
        /// This causes same DB assigned on multiple nodes.
        /// This causes same DB assigned on multiple nodes.
        ///
        nodeIdx = 0;
        cnt = 0;
        for (int i = 0; i < numNodeRemain; ++i) {
            //mDbCoreNum[vDbChunkNames[i]] = nodeIdx;
            mDbCoreNum.insert(pair <string, int> (vDbChunkNames[i], nodeIdx)); 
            if (cnt == numDBPerNode - 1) {
                nodeIdx++;
                cnt = 0;
            }
            else cnt++;
        }
        
        /// DEBUG
        cout << "mDbCoreNum.size() = " << (int) mDbCoreNum.size() << endl;
        for (std::map<std::string, int>::iterator im = mDbCoreNum.begin(); im != mDbCoreNum.end(); ++im) {
            std::cout << "mDbCoreNum = " << im->first << " " << im->second << std::endl;
        }
        ///
        
        gf.mDbCoreNum = &mDbCoreNum;
        vDbChunkNames.clear();
        */
    }
    MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////
    ///
    /// Now vWorkItems has all the work items in the form of 
    /// (query_file_name,DB_chunk_name)
    ///
    
    uint64_t numFilesToCreate = 0;
    uint64_t numWorkItemsForEach = 0;
    if (atoi(argv[3]) > 0) {
        if (MPI_myId == 0) {
            ///
            /// Divide the whole work items into several work item lists (stored
            /// in .
            /// Each iteration consists of calling blast and saving output
            /// into an output file.
            ///
            //uint64_t numQueryBunch = str2uint(string(argv[3]));
            uint64_t numQueryBunch = atoi(argv[3]);
            if (numQueryBunch > numQueryFiles) {
                cerr << "### ERROR: should be numQueryBunch <= numQueryFiles\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            //if (numQueryFiles % numQueryBunch) {
                //cerr << "### ERROR: should be numQueryFiles % numQueryBunch == 0\n";
                //MPI_Abort(MPI_COMM_WORLD, 1);
            //}
            numFilesToCreate = numQueryFiles / numQueryBunch;
            uint64_t numTotalWorkItems = vWorkItems.size();
            numWorkItemsForEach = numQueryBunch * numDbChunks;
            uint64_t numRemains = numTotalWorkItems % numWorkItemsForEach;
            //if (numQueryBunch % numDbChunks) {
                //cerr << "### ERROR: should be numQueryBunch % numQueryFiles == 0\n";
                //MPI_Abort(MPI_COMM_WORLD, 1);
            //}
            
            uint64_t k = 0;
            uint64_t i = 0;
            for (; i < numFilesToCreate; ++i) {
                string newFileName = prefix + "-workitems-" + uint2str(i) + ".txt";
                ofstream workItemFile(newFileName.c_str());
                cout << "### INFO: Work item file name (" << i << ") = " 
                     << newFileName << endl;
                if (workItemFile.is_open()) {
                    for (uint64_t j = 0; j < numWorkItemsForEach; ++j, ++k) {
                        workItemFile << vWorkItems[k] << endl;
                    }
                    workItemFile.close();
                }
                else {
                    cerr << "### ERROR: work item file open error.\n";
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            if (numRemains) {
                string newFileName = prefix + "-workitems-" + uint2str(i) + ".txt";
                ofstream workItemFile(newFileName.c_str());
                cout << "### INFO: Work item file name (remains) = " 
                     << newFileName << endl;
                if (workItemFile.is_open()) {
                    for (uint64_t j = 0; j < numRemains; ++j, ++k) {
                        workItemFile << vWorkItems[k] << endl;
                    }
                    workItemFile.close();
                    numFilesToCreate++; /// for one remains file
                }
                else {
                    cerr << "### ERROR: work item file open error.\n";
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            cout << "### INFO: num work item files = " << numFilesToCreate << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////


        ///
        /// Iteratively call blast and save results for numFilesToCreate
        /// times.
        ///

        /// Broadcast the number of sub work item files
        MPI_Bcast(&numFilesToCreate, 1, MPI_INT, 0, MPI_COMM_WORLD);

        //int x = 0;
        //int y = 0; 
        for (uint64_t n = 0; n < numFilesToCreate; ++n) {
            vector<string> vSubWorkitems;
            string newWorkItemFileName = prefix + "-workitems-" + uint2str(n) + ".txt";
            /*
            if (MPI_myId == 0) {       
                y = x + numWorkItemsForEach;         
                cout << "### INFO: Start processing " << newWorkItemFileName
                     << endl;
                if (vWorkItems.begin()+x+y <= vWorkItems.end())
                    copy(vWorkItems.begin()+x, vWorkItems.begin()+x+y, vSubWorkitems.begin());
                else 
                    copy(vWorkItems.begin()+x, vWorkItems.end(), vSubWorkitems.begin());
                cout << "vSubWorkitems.size() = " << vSubWorkitems.size() << endl;
                x += numWorkItemsForEach;
            }
            */
            
            /// DEBUG
            //char *temp = "workitems0.txt";
            //nvecRes = mr2->map(temp, &run_blast2, &gf);
            ///

            /// ////////////////////////////////////////////////////////
            #ifdef DEBUG
            c2 = clock();
            t2 = time(NULL);
            #endif
            nvecRes = mr2->map((char*)newWorkItemFileName.c_str(), &run_blast2, &gf);
            //nvecRes = mr2->map(vSubWorkitems, &run_blast2, &gf);
            #ifdef DEBUG
            c3 = clock();
            t3 = time(NULL);
            accMapWTime += (long)t3-t2;
            accMapCTime += (long)c3-c2;
            #endif
            
            if (MPI_myId == 0) {
                cout << "### INFO: [Rank " << MPI_myId << "] Num res from map() = "
                     << nvecRes << endl;
            }
            
            ///
            /// Old version
            ///
            //mr2->collate(NULL);
            //uint64_t nunique = mr2->reduce(&remove_DBname_from_key, NULL);
            //MPI_Barrier(MPI_COMM_WORLD);
            //mr2->gather(1);
            ////mr2->sort_keys(&key_compare_int); /// Really need this?
            //mr2->convert();
            //nunique = mr2->reduce(&sort_multivalues_by_score, NULL);
            
            /// 
            /// New version with sorting done in workers
            /// 
            mr2->collate(NULL);
            uint64_t nunique = mr2->reduce(&remove_DBname_from_key, NULL);
            mr2->collate(NULL);
            nunique = mr2->reduce(&sort_multivalues_by_score2, NULL);
            
            ///
            /// 10.20.2010
            /// Try to save output file at each worker. 
            /// Note: No need to call gather() now.
            ///       The rank-0 should print out its KV in a file.
            /// 
            //if (MPI_myId != 0) {
                string workerOutputFileName = prefix + "-" + uint2str(n) + "-" + uint2str(MPI_myId) + ".txt";
                FILE *outFile = fopen(workerOutputFileName.c_str(), "w");
                if (outFile) mr2->kv->print2file2(1, 5, 5, outFile);
                else {
                    cerr << "### ERROR: outFile open error.\n";
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                fclose(outFile);
                cout << "### INFO: File saved! => " << workerOutputFileName << endl; 
            //}
            
            ///
            /// Save history
            ///
            if (MPI_myId == 0) {                
                //char *histFileName = "history.txt";
                string histFileName = prefix + "-history.txt";
                FILE *histFile = fopen(histFileName.c_str(), "a");
                time_t timer;
                timer = time(NULL);
                fprintf(histFile, "%s,%s", newWorkItemFileName.c_str(), asctime(localtime(&timer)));
                fclose(histFile);
            }
                
            MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////            
            
            /*
            mr2->gather(1);
            /// ////////////////////////////////////////////////////////
            
            if (MPI_myId == 0) {
                //cout << "### INFO: Saving output...\n";
                string outFileName = prefix  + "-" + uint2str(n) + ".txt";
                //cout << "### INFO: outFileName = " << outFileName << endl;
                FILE *outFile = NULL;
                //if (n == 0) outFile = fopen(outFileName.c_str(), "w");
                //else outFile = fopen(outFileName.c_str(), "a");
                outFile = fopen(outFileName.c_str(), "w");
                if (outFile) mr2->kv->print2file2(1, 5, 5, outFile);
                else { 
                    cerr << "### ERROR: outFile open error.\n";
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                cout << "### INFO: File saved! => " << outFileName << endl;  
                
                /// save history
                //char *histFileName = "history.txt";
                string histFileName = prefix + "-history.txt";
                FILE *histFile = fopen(histFileName.c_str(), "a");
                time_t timer;
                timer = time(NULL);
                fprintf(histFile, "%s,%s\n", newWorkItemFileName.c_str(), asctime(localtime(&timer)));
                fclose(histFile);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            */
        }
    }
    else { /// Process whole work items
        //char *workItemFileName = prefix + "-workitems.txt";
        string workItemFileName = prefix + "-workitems.txt";
        ///
        /// Saving work items in a file just for history.
        ///
        if (MPI_myId == 0) {
            ofstream workItemFile(workItemFileName.c_str());
            if (workItemFile.is_open()) {
                for (uint64_t i = 0; i < vWorkItems.size(); ++i)
                    workItemFile << vWorkItems[i] << endl;
                workItemFile.close();
            }
            else {
                cerr << "### ERROR: workitems.txt not found.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////

        /// DEBUG
        //char* workItemFileName2 = "test_workitems.txt";
        //nvecRes = mr2->map(workItemFileName2, &run_blast2, &gf);
        ///

        ///
        /// MAPREDUCE PROCEDURE
        ///
        /// 1. map(): KV, each rank get its work item in "query_file,DB_id" form
        /// 2. collate()/convert(): KV->KMV. NOTE: collate() can be called. But no
        ///    need to aggregate the same KV from other ranks.
        /// 3. reduce(): KMV->KV. Remove DB_ids from keys.
        /// 4. gather(): KV. Collect all KVs across all ranks.
        /// 5. sort_keys(): KV. Sort keys in the collected KV.
        /// 6. convert(): KV->KMV. All blast results for each query_id is collected.
        /// 7. reduce(): KMV->KV. Sort multivalues by hit score.
        /// 8. mr->kv->print2file(): Print output.
        ///
         
        ///
        /// New version with sorting done in workders
        /// 1. map(): blast call with assigned work item
        /// 2. collate(): creates KMV with (Qid,DBid) as key
        /// 3. reduce(): convert <(Qid,DBid), {blast_result(s)}> -> <Qid, {blast_result(s)}>
        /// 4. collate(): creates KMV with Qid as key
        /// 5. reduce(): sort multivalues by bit score
        /// 6. gather(): gather KVs to master
        /// 7. print(): printing
        ///

        /// 1.
        #ifdef DEBUG
        c2 = clock();
        t2 = time(NULL);
        #endif
        /// Process work items in the file, workitems.txt
        //nvecRes = mr2->map(workItemFileName, &run_blast2, &gf);
        /// Process work items in vWorkItems
        nvecRes = mr2->map(vWorkItems, &run_blast2, &gf);
        #ifdef DEBUG
        c3 = clock();
        t3 = time(NULL);
        accMapCTime += (long)c3-c2;
        accMapWTime += (long)t3-t2;
        #endif
        
        if (MPI_myId == 0) {
            cout << "### [Rank " << MPI_myId << "] Num results from map() = "
                 << nvecRes << endl;
        }
        //mr2->print(-1, 1, 5, 5);

        /// 2.
        mr2->collate(NULL);
        //mr2->convert();
        //mr2->print(-1, 1, 5, 5);

        /// 3.
        uint64_t nunique = mr2->reduce(&remove_DBname_from_key, NULL);
        //mr2->print(-1, 1, 5, 5);

        /// 3-1.
        mr2->collate(NULL); 
        //mr2->print(-1, 1, 5, 5);
        
        //nunique = mr2->reduce(&remove_DBname_from_key, NULL);
        nunique = mr2->reduce(&sort_multivalues_by_score2, NULL);
        //mr2->print(-1, 1, 5, 5);
        
        ///
        /// 10.20.2010
        /// Try to save output file at each worker. 
        /// Note: No need to call gather() now.
        ///       The rank-0 should print out its KV in a file.
        /// 
        //if (MPI_myId != 0) {
            string workerOutputFileName = prefix + "-" + uint2str(MPI_myId) + ".txt";
            //cout << "### INFO: outFileName = " << outFileName << endl;
            FILE *outFile = fopen(workerOutputFileName.c_str(), "w");
            if (outFile) mr2->kv->print2file2(1, 5, 5, outFile);
            else {
                cerr << "### ERROR: outFile open error.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            fclose(outFile);
            cout << "### INFO: File saved! => " << workerOutputFileName << endl; 
        //}
        
        MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////////
        /*
        /// 4.
        mr2->gather(1);
        //mr2->print(0, 1, 5, 5);

        /// 5. 
        /// Note: if you want to keep the original query order.
        //mr2->sort_keys(&key_compare_int); /// Really need this?
        //mr2->print(0, 1, 5, 5);

        /// 6.
        //mr2->convert();
        //mr2->print(0, 1, 5, 5);

        /// 7.
        /// Note: sort_multivalues() cannot be used in this implementation
        /// because KMV could have a form
        /// of {"str str", str} for Blast result strings.
        ///
        //nunique = mr2->reduce(&sort_multivalues_by_score, NULL);  
        //mr2->print(0, 1, 5, 5);

        /// 8. save output
        if (MPI_myId == 0) {
            //cout << "Saving...\n";
            string outFileName = string(argv[2]) + ".txt";
            //cout << "### INFO: outFileName = " << outFileName << endl;
            FILE *outFile = fopen(outFileName.c_str(), "w");
            if (outFile) mr2->kv->print2file2(1, 5, 5, outFile);
            else {
                cerr << "### ERROR: outFile open error.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            fclose(outFile);
            cout << "### INFO: File saved! => " << outFileName << endl;  
        }
        //string outFileName = string(argv[2]) + ".txt";
        //cout << "### INFO: outFileName = " << outFileName << endl;
        ////FILE *outFile = fopen(outFileName.c_str(), "w");
        //gf.outFileName = (char*)outFileName.c_str();
        //mr2->map(mr2, &save_output, &gf);   
        //cout << "### INFO: File saved!\n";   
        */
    }
    delete mr2;
    
//} catch (MPI::Exception e) {
    //std::cout << "MPI ERROR: " << e.Get_error_code() 
              //<< " -" << e.Get_error_string() 
              //<< std::endl;
    //MPI_Abort(MPI_COMM_WORLD, 1);
//}              

    MPI_Finalize();
    
    #ifdef DEBUG
    t1 = time(NULL);
    c1 = clock();
    //printf ("### TIME: [Rank %d] elapsed total wall clock time,%ld\n", MPI_myId, (long) (t1 - t0));
    //printf ("### TIME: [Rank %d] elapsed map wall clock time,%ld\n", MPI_myId, (long) (accMapWTime));
    //printf ("### TIME: [Rank %d] elapsed total CPU time,%f\n", MPI_myId, (double) (c1 - c0)/CLOCKS_PER_SEC);
    //printf ("### TIME: [Rank %d] elapsed map CPU time,%f\n", MPI_myId, (double) (accMapCTime)/CLOCKS_PER_SEC);    
    
    /// FORMAT: rank,totalRuntime,mapRuntime,totalClock,mapClock
    printf ("### TIME: rank,totalRuntime,mapRuntime,totalClock,mapClock,%d,%ld,%ld,%.2f,%.2f\n", MPI_myId, 
            (long) (t1 - t0),
            (long) (accMapWTime),
            (double) (c1 - c0) / CLOCKS_PER_SEC,
            (double) (accMapCTime) / CLOCKS_PER_SEC);
    #endif
    
    return 0;
}

/** MR-MPI Map function - settup NCBI C++ Toolkit env and call blast
 * @param itask
 * @param file
 * @param kv
 * @param ptr
 */
void run_blast2(int itask, char *file, KeyValue *kv, void *ptr)
{
    #ifdef DEBUG
    time_t t0, t1, t2, t3, accBlastCallTime=0;
    time_t t4, t5, accDBLoadingTime=0;
    t0 = time(NULL);
    #endif 
    
    giftBox *gf = (giftBox *) ptr;
    
    /// 
    /// Make a option handle and cblastinputsource
    ///
    EProgram program = ProgramNameToEnum("blastn");
    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));
    set_default_opts(opts);
    opts->Validate();

    /// DEBUG
    if (gf->myId == 0)
        opts->GetOptions().DebugDumpText(cerr, "opts", 1);
    ///

    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
        throw std::runtime_error("Could not initialize object manager");
    }

    bool isProtein = false;
    SDataLoaderConfig dlconfig(isProtein);
    CBlastInputSourceConfig iconfig(dlconfig);
    //fprintf(stdout, "### INFO: [Rank %d]: %s, itask = %d, Work item = %s\n",
    //gf->myId, gf->pName, itask, file);
    //string qFileName = "test2.query";
    //ifstream queryFile(qFileName.c_str());
    //CBlastFastaInputSource fasta_input(qFileName, iconfig);

    ///
    /// Split work item => query file + db chunk name
    ///
    string workItem(file);
    //gf->workItemProcessed = workItem;
    vector<string> vWorkItem = split(workItem, ',');
    ifstream queryFile(vWorkItem[0].c_str());
    if (!queryFile.is_open()) {
        cerr << "### ERROR: queryFile open error.\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    //char *key_tokens = strtok(key, ',');

    ///
    /// Use query string instead of file input
    ///
    string header, seq;
    vector<string> vecSeq;

    ///
    /// Target db chunk setting
    ///
    string dbChunkName = vWorkItem[1];
    #ifdef DEBUG
    t4 = time(NULL);
    #endif
    const CSearchDatabase target_db(dbChunkName,
                                    CSearchDatabase::eBlastDbIsNucleotide);
    
    #ifdef DEBUG
    t5 = time(NULL);
    accDBLoadingTime += (long)t5-t4;
    #endif
    
    ///
    /// Read seq(s) and run blast
    ///
    while (!getline(queryFile, header).eof()) {
        getline(queryFile, seq);

        vector<string> vHeaders;

        if (header.length() == 0 || seq.length() == 0) {
            fprintf(stderr, "### ERROR: [Rank %d]: %s, itask = %d, file = %s, \
                    Seq read error.\n", gf->myId, gf->pName, itask, file);
            MPI_Abort(MPI_COMM_WORLD, 1);       
            exit(0);
        }
        //fprintf(stdout, "### INFO: [Rank %d]: %s, itask = %d, %s\n",
        //gf->myId, gf->pName, itask, header.c_str());
        
        /// NOTE: This should be removed!
        string query = header + '\n' + seq;
        vHeaders.push_back(header);
        ///

        ///
        /// ACCUMMULATE SEQS WHEN MAXQUERYNUM > 1
        ///
        uint64_t j = 1;
        while (j < MAXQUERYNUM && !getline(queryFile, header).eof()) {
            getline(queryFile, seq);
            //fprintf(stdout, "### INFO: [Rank %d]: %s, itask = %d, %s, \
            //num queries = %d\n",
            //gf->myId, gf->pName, itask, header.c_str(), j + 1);
            query = query + '\n' + header + '\n' + seq;
            vHeaders.push_back(header);
            j++;
        }

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
        #ifdef DEBUG
        t2 = time(NULL);
        #endif
        CSearchResultSet results = *blaster.Run();
        #ifdef DEBUG
        t3 = time(NULL);
        accBlastCallTime += (long)t3-t2;
        #endif
        
        ///
        /// GET WARNING MESSAGES
        ///
        for (uint64_t i = 0; i < results.GetNumResults(); i++) {
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

        ///
        /// GET THE RESULTS
        ///
        //cout << "results.GetNumResults() = " << results.GetNumResults() << endl;
        fprintf(stdout, "### INFO: [Rank %d]: %s, itask = %d, work item = %s, \
                Num query = %d, Num Res = %d\n",
                gf->myId, gf->pName, itask, file, results.GetNumQueries(),
                results.GetNumResults());
        //cout << "results.GetNumQueries() = " << results.GetNumQueries() << endl;

        //string exlusionFileName = "excluded-rank"+int2str(gf->myId)+".txt";
        //ofstream exclusionSaveFile(exlusionFileName.c_str());
        
        for (uint64_t i = 0; i < results.GetNumResults(); i++) {

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

                    uint64_t alignLen, qStart, qEnd, sStart, sEnd;
                    //CRange<TSeqPos> qRange, sRange;
                    //qRange    = s.GetSeqRange(QUERY);
                    //sRange    = s.GetSeqRange(SUBJECT);
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
                    //uint64_t origLen       = str2uint(vQeuryId[3]);   /// LENGTH OF THE ORIG SEQ
                    int qCutLocStart         = str2int(vQeuryId[4]);    /// CUT COORDINATES - START
                    uint64_t qCutLocEnd      = str2uint(vQeuryId[5]);   /// CUT COORDINATES - END
                    //fprintf(stderr, "### INFO: [Rank %d]: %s, itask = %d, %s %s %d %d %d\n",
                    //gf->myId, gf->pName, itask, qGi.c_str(), qid.c_str(),
                    //origLen, qCutLocStart, qCutLocEnd);

                    /// 
                    /// ADD A CSV BLAST RESULT TO KV
                    /// 
                    if (!check_exclusion(qGi, seqId, qCutLocStart, qCutLocEnd, 
                        sStart, sEnd, EXCLUSION)) {
                        char blastRes[MAXSTR];

                        /// DEBUG
                        //fprintf(stdout, "%s,%s,%s,%d,%d,%d,%1.e,%d\n",
                        //qid.c_str(), qHeader.substr(1, MAXSTR).c_str(), seqId.c_str(),
                        //alignLen, qStart, qEnd,
                        ////eValue, genericScore, bitScore);
                        //eValue, bitScore);
                        ///

                        //sprintf(blastRes, "%s,%s,%d,%d,%d,%1.e,%d",
                                //qHeader.substr(0, MAXSTR).c_str(), seqId.c_str(),
                                //alignLen, qStart, qEnd,
                                ////eValue, genericScore, bitScore);
                                //eValue, bitScore);
                        
                        /// 
                        /// 10.9.2010
                        /// Add db chunk name and query file name at the end of 
                        /// blast result string for recording history.
                        /// 
                        //sprintf(blastRes, "%s,%s,%d,%d,%d,%1.e,%d,%s",
                                //qHeader.substr(0, MAXSTR).c_str(), seqId.c_str(),
                                //alignLen, qStart, qEnd,
                                ////eValue, genericScore, bitScore);
                                //eValue, bitScore, file);
                         sprintf(blastRes, "%s,%s,%d,%d,%d,%1.e,%d",
                                qHeader.substr(0, MAXSTR).c_str(), seqId.c_str(),
                                alignLen, qStart, qEnd,
                                eValue, bitScore);

                        char *newKey = (char*)((qid + ',' + dbChunkName).c_str());

                        /// 
                        /// ADD <KEY = "QUERYID,CHUNKNAME", VALUE="BLASTvecResULT">
                        /// TO KV
                        /// 
                        kv->add(newKey, strlen(newKey) + 1, blastRes, strlen(blastRes) + 1);
                        
                        /// 
                        /// 10.8.2010
                        /// The query file name is added to newKey for recording 
                        /// history.
                        /// 
                        //vector<string> vQueryFileNameTokens = split(vWorkItem[0], '.');
                        //newkey = 
                        //kv->add(newKey, strlen(newKey) + 1, blastRes, strlen(blastRes) + 1);
                    }
                    else {
                        fprintf(stdout, "### INFO: EXCLUSION - qGi=%s sGi=%s \
                                qCutLocStart=%d, qCutLocEnd=%d, sStart=%d, sEnd=%d\n",
                                qGi.c_str(), seqId.c_str(), qCutLocStart, qCutLocEnd, 
                                sStart, sEnd);

                        //if (exclusionSaveFile.is_open()) {
                        //char excluded[MAXSTR];
                        //sprintf(stderr, "qGi=%s sGi=%s \
                        //qCutLocStart=%d, qCutLocEnd=%d, sStart=%d, sEnd=%d\n",
                        //qGi.c_str(), seqId.c_str(), qCutLocStart, qCutLocEnd, sStart, sEnd);
                    }
                }
            }
        }
    } /// END WHILE
    
    #ifdef DEBUG
    t1 = time(NULL);
    /// FORMAT: rank,total_runblast2,blaster.Run,DBLoadingTime
    printf ("### TIME (in runblast2): rank,total_runblast2,blaster.Run,DBLoadingTime,%d,%ld,%ld,%ld\n", 
            gf->myId, (long)(t1 - t0), (long)accBlastCallTime, (long)accDBLoadingTime);
    #endif
}

///* ---------------------------------------------------------------------- */
//void save_output(uint64_t itask, char *key, int keybytes, char *value,
        //int valuebytes, KeyValue *kv, void *ptr)
///* ---------------------------------------------------------------------- */
//{
    //vector<string> vResults = split(string(value), '>');
    //cout << "vResults = " << vResults.size() << endl;
    
    //giftBox *gf = (giftBox *) ptr;
    //FILE *outFile = fopen(gf->outFileName, "w");
    
    //for (uint64_t i = 0; i < vResults.size(); ++i) 
        //fprintf(outFile, "%s\n", (char*)vResults[i].c_str());
    
    //fclose(outFile);
//}

/** Sort function - Passed to MR-MPI sort_values() for sorting blast result
 * string by bit score.
 * @param str1
 * @param str2
 */
bool mysort(string str1, string str2)
{
    vector<string> vStr1Tokens = split(str1, ',');
    vector<string> vStr2Tokens = split(str2, ',');
    //assert(vStr1Tokens.size() == 8);
    //assert(vStr2Tokens.size() == 8);
    //cout << "vStr1Tokens.size() = " << vStr1Tokens.size() << endl;
    //cout << "vStr2Tokens.size() = " << vStr2Tokens.size() << endl;
    /// sort by bit score
    //cout << "vStr1Tokens[6] = " << vStr1Tokens[6] << endl;
    //cout << "vStr2Tokens[6] = " << vStr2Tokens[6] << endl;
    return (str2uint(vStr1Tokens[6]) > str2uint(vStr2Tokens[6]));
}

/** Sort by bit score - Passed to MR-MPI reduce() for sorting KMVs by bit score.
 * @param key
 * @param keybytes
 * @param multivalue: collected blast result strings. There could be more than
 * nvalues results in it. The separator '>' should be used for splitting the 
 * resutls.
 * @param nvalues
 * @param valuebytes
 * @param kv
 * @param ptr
 */
void sort_multivalues_by_score2(char *key, int keybytes, char *multivalue,
                   int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
    /// 
    /// Note: Even if nvalues == 1, multivalue could have multiple
    /// blast results.
    /// 
    uint64_t sum = 0;
    uint64_t idx = 0;
    string newMultiValue;
    for (uint64_t i = 0; i < nvalues; ++i) {
        //if (i > 0) newMultiValue += " ";
        sum += valuebytes[i];
        char temp[valuebytes[i]];
        for (uint64_t j = 0; j < valuebytes[i]; ++j) {
            temp[j] = multivalue[idx];
            idx++;
        }
        newMultiValue += string(temp);
    }

    /// 
    /// SORTING
    /// 
    vector<string> vTokens = split(newMultiValue, '>');
    //cout << "vTokems from newMultiValue = " << vTokens.size() << endl;
    if (vTokens.size() > 2) {
        //for (int i = 0; i < vTokens.size(); ++i) {
        //cout << "vTokens[i] = " << vTokens[i] << endl;
        //}
        vTokens.erase(vTokens.begin()); /// There always is a null token.
        sort(vTokens.begin(), vTokens.end(), mysort);
    }
    else if (vTokens.size() == 2) {     /// If there is on one line of result.
        vTokens.erase(vTokens.begin()); /// Just erase the first blank token from vTokens.
    }
    
    string newMultiValue2;
    for (int i = 0; i < vTokens.size(); ++i) {
        ////cout << vTokens[i] << endl;
        //kv->add(key, strlen(key) + 1, (char*)vTokens[i].c_str(), 
                //vTokens[i].size() + 1);
        newMultiValue2 += '>' + vTokens[i];
    }
    kv->add(key, strlen(key) + 1, (char*)newMultiValue2.c_str(), newMultiValue2.size() + 1);
}

/** Sort by bit score - Passed to MR-MPI reduce() for sorting KMVs by bit score.
 * @param key
 * @param keybytes
 * @param multivalue: collected blast result strings. There could be more than
 * nvalues results in it. The separator '>' should be used for splitting the 
 * resutls.
 * @param nvalues
 * @param valuebytes
 * @param kv
 * @param ptr
 */
/*
void sort_multivalues_by_score(char *key, int keybytes, char *multivalue,
                   int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
    /// Note: Even if nvalues == 1, multivalue could have multiple
    /// blast results.
    uint64_t sum = 0;
    uint64_t idx = 0;
    string newMultiValue;
    for (uint64_t i = 0; i < nvalues; ++i) {
        //if (i > 0) newMultiValue += " ";
        sum += valuebytes[i];
        char temp[valuebytes[i]];
        for (uint64_t j = 0; j < valuebytes[i]; ++j) {
            temp[j] = multivalue[idx];
            idx++;
        }
        newMultiValue += string(temp);
    }

    /// Sort
    vector<string> vTokens = split(newMultiValue, '>');
    //cout << "vTokems from newMultiValue = " << vTokens.size() << endl;
    if (vTokens.size() > 2) {
        //for (int i = 0; i < vTokens.size(); ++i) {
        //cout << "vTokens[i] = " << vTokens[i] << endl;
        //}
        vTokens.erase(vTokens.begin()); /// There always is a null token.
        sort(vTokens.begin(), vTokens.end(), mysort);
    }
    else if (vTokens.size() == 2) {     /// If there is on one line of result.
        vTokens.erase(vTokens.begin()); /// Just erase the first blank token from vTokens.
    }

    for (int i = 0; i < vTokens.size(); ++i) {
        //cout << vTokens[i] << endl;
        kv->add(key, strlen(key) + 1, (char*)vTokens[i].c_str(), 
                vTokens[i].size() + 1);
    }
}
*/

/** A Reduce user-defined function - Passed to MR-MPI reduce() for converting
 * (Qid,DBid) key into Qid for further aggregating.
 * @param key
 * @param keybytes
 * @param multivalue: collected blast result strings.  
 * @param nvalues
 * @param valuebytes
 * @param kv
 * @param ptr
 */
void remove_DBname_from_key(char *key, int keybytes, char *multivalue,
                     int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
    uint64_t sum = 0;
    uint64_t idx = 0;
    string newMultiValue;

    /// 
    /// Note: passing multivalue in the function input argument to
    /// kv->add is not working.
    /// Force to copy multivalue to temp and pass temp to kv->add
    /// 
    for (uint64_t i = 0; i < nvalues; ++i) {
        if (i > 0) newMultiValue += " ";
        sum += valuebytes[i];
        char temp[valuebytes[i]];
        for (uint64_t j = 0; j < valuebytes[i]; ++j) {
            temp[j] = multivalue[idx];
            idx++;
        }
        newMultiValue += string(temp);
    }
    char *newKey = strtok(key, ",");
    kv->add(newKey, strlen(newKey) + 1, (char*)newMultiValue.c_str(), sum);
}

/** Check exclusion - Based on the coordinates of query and subject, decide 
 * whether the result should be included in the final result or not.
 * @param qGi: query ID
 * @param sGi: subject ID
 * @param qCutLocStart: Cutting start location of the query from original 
 * refseq_genomic fasta input.
 * @param qCutLocEnd: Cutting end location of the query from original 
 * refseq_genomic fasta input.
 * @param sStart: subject alignment start location.
 * @param sEnd: subject alignemnt end location.
 * @param threshold: overlap threshold (default = 100bp).
 */
bool check_exclusion(string qGi, string sGi, int qCutLocStart,
                     uint64_t qCutLocEnd, uint64_t sStart, uint64_t sEnd,
                     uint64_t threshold)
{
    /// 
    /// To exclude Blast result from the original sequence from which
    /// the input query is originated (sampled). Basically if qGi == sGi
    /// and qCutLocStart is similar with sStart and qCutLocEnd is similar
    /// with sEnd in terms of coordinates, the result should be excluded.
    /// 
    /// Orig seq: ----------------XXXXXXXXXXXXXXXX----------------------
    ///                           |              |
    ///                     qCutLocStart      qCutLocEnd
    /// 
    /// Query:                    XXXXXXXXXXXXXXXX
    ///                              |          |
    ///                           qStart       qEnd
    /// 
    /// Subject:   ------------------XXXXXXXXXXXX-----------------------
    ///                              |          |
    ///                           sStart      sEnd
    /// 
    bool ret = false;

    if (qGi == sGi) {
        if (qCutLocStart < 0) {
            /// 
            /// In >gi|222299657|18|3605|-400|3604
            /// -400|3604 means query[-400:3604] in Python.
            /// 
            qCutLocStart = qCutLocEnd + 1 - qCutLocStart;
            qCutLocEnd += 1;
        }
        if ((qCutLocStart - threshold <= sStart 
            && sStart <= qCutLocStart + threshold) 
            && (qCutLocEnd - threshold <= sEnd 
            && sEnd <= qCutLocEnd + threshold))
            ret = true;
    }

    return ret;
}

/** Set Blast options
 * @param optsHandle
 */
void set_default_opts(CRef<CBlastOptionsHandle> optsHandle)
{
    optsHandle->SetEvalueThreshold(EVALUE);
    //optsHandle->SetMatchReward(0);
    //optsHandle->SetMismatchPenalty(0);
    //optsHandle->SetMatrixName("BLOSUM62");
    optsHandle->SetCutoffScore(CUTOFFSCORE);
    if (CBlastNucleotideOptionsHandle* nuclHandle =
                dynamic_cast<CBlastNucleotideOptionsHandle*>(&*optsHandle)) {

        //nuclHandle->SetMatchReward(0);
        //nuclHandle->SetMismatchPenalty(0);
        nuclHandle->SetMatrixName("BLOSUM62");
        nuclHandle->SetCutoffScore(CUTOFFSCORE);
        nuclHandle->SetEvalueThreshold(EVALUE);
    }

    return;
}

/*
int key_compare_str(char *p1, int len1, char *p2, int len2)
{
    return strcmp(p1, p2);
}
*/

/*
int key_compare_str2(char *p1, int len1, char *p2, int len2)
{
    ///  
    /// Key = (query_id,DB_id)
    /// Fist, compare query_ids. If those are same, compare
    /// DB_ids.
    /// 
    string str1(p1);
    string str2(p2);
    vector<string> vStr1 = split(str1, ',');
    vector<string> vStr2 = split(str2, ',');
    uint64_t i1 = str2uint(vStr1[0]);
    uint64_t i2 = str2uint(vStr2[0]);
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
*/

/*
int key_compare_int(char *p1, int len1, char *p2, int len2)
{
    //int i1 = *(int *) p1;
    //int i2 = *(int *) p2;
    uint64_t i1 = str2uint(string(p1));
    uint64_t i2 = str2uint(string(p2));
    /// Ascending order
    if (i1 > i2) return 1;
    else if (i1 < i2) return -1;
    /// Descending order
    //if (i1 > i2) return -1;
    //else if (i1 < i2) return 1;
    else return 0;
}
*/


/*
void get_node_name(int itask, KeyValue *kv, void *ptr)
{
    procName *nn = (procName *) ptr;
    char id[MAXSTR];
    sprintf(id, "%d", nn->myId);
    kv->add(nn->pName, strlen(nn->pName) + 1, id, strlen(id) + 1);
}

void collect_node_names(char *key, int keybytes, char *multivalue,
                        int nvalues, int *valuebytes, KeyValue *kv,
                        void *ptr)
{
    kv->add(key, keybytes, multivalue, sizeof(multivalue));
}

void save_node_names(uint64_t itask, char *key, int keybytes, char *value,
                     int valuebytes, KeyValue *kv, void *ptr)
{
    //multimap<string, uint64_t> *nNames = (multimap<string, uint64_t> *) ptr;
    vector<string>* vecNames = (vector<string> *) ptr;
    string pName(key);
    //uint64_t myid = atoi(value);
    //nNames->insert(pair<string, uint64_t>(pName, myid));
    vector<string>::const_iterator loc = find(vecNames->begin(), vecNames->end(), 
                                                pName);
    if (loc == vecNames->end())
        vecNames->push_back(pName);
}
*/

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
/*
string int2str(int number)
{
    stringstream ss;
    ss << number;
    return ss.str();
}
* */

string uint2str(uint64_t number)
{
    stringstream ss;
    ss << number;
    return ss.str();
}

uint64_t str2uint(string str)
{
    std::stringstream ss;
    ss << str;
    uint64_t f;
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


/* 
void run_blast(int itask, char *file, KeyValue *kv, void *ptr)
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
            fprintf(stderr, "### ERROR: [Rank %d]: %s, itask = %d, file = %s, Seq read error.\n",
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
        uint64_t j = 1;
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
        //uint64_t oLen              = str2uint(vQeuryId[3]);   /// LENGTH OF THE ORIG SEQ
        //int qStart               = str2int(vQeuryId[4]);    /// CUT COORDINATES - START
        //uint64_t qEnd              = str2uint(vQeuryId[5]);   /// CUT COORDINATES - END
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
        for (uint64_t i = 0; i < results.GetNumResults(); i++) {
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
        for (uint64_t i = 0; i < results.GetNumResults(); i++) {

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

                    uint64_t alignLen, qStart, qEnd, sStart, sEnd;
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
                    uint64_t origLen           = str2uint(vQeuryId[3]);   /// LENGTH OF THE ORIG SEQ
                    int qCutLocStart         = str2int(vQeuryId[4]);    /// CUT COORDINATES - START
                    uint64_t qCutLocEnd        = str2uint(vQeuryId[5]);   /// CUT COORDINATES - END
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
                        fprintf(stderr, "### INFO: EXCLUSION - qGi=%s sGi=%s qCutLocStart=%d, qCutLocEnd=%d, sStart=%d, sEnd=%d\n",
                            qGi.c_str(), seqId.c_str(), qCutLocStart, qCutLocEnd, sStart, sEnd);

                }
            }
        }
    } // END WHILE
}
*/

/** atexit() user-defined function - to ensure MPI_Finalize()
 */
 /*
void AtExit(void)
{
    if (!MPI_Finalize()) MPI_Abort(MPI_COMM_WORLD, 1);
}
*/

/// EOF

