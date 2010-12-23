////////////////////////////////////////////////////////////////////////
//
//  MPI-MR-BLAST: Parallelizing BLAST on MR-MPI
//
//  Author: Seung-Jin Sul
//          (ssul@jcvi.org)
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
#include "mrmpi/mapreduce.h"
#include "mrmpi/keyvalue.h"

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
#include <string.h>

///For CSeq_align for output
#include <objects/seqalign/Seq_align.hpp>
#include <util/range.hpp>

/// For typedef unsigned long long int uint64_t
#include <stdint.h>

/// For timing
#include <sys/time.h>
//#define DEBUG 1

/// For queue
#include <queue>

/// Conf file processing from mrblast.conf
#include "conf/ConfigFile.h"

USING_NCBI_SCOPE;
USING_SCOPE(blast);




const uint64_t MAXSTR = 256;
unsigned int  MAXQUERYNUM = 0; /// num query to accumulate from qeury file
const uint64_t QUERY = 0;
const uint64_t SUBJECT = 1;
//const uint64_t EXCLUSION = 100;     /// Exclusion threshold = 100bp

/// Blast options
unsigned int CUTOFFSCORE = 0;   
double EVALUE = 0.0;

string DBCHUNKLISTFILE = "dbchunks.txt";
string BLASTOPTS = "";
unsigned int NUMTOTALDBCHUNKS = 0;
int ITER = 0;

/// Blast DB setting
static CSearchDatabase *pTargetDb = 0;
string prevDbChunkName;

/// Misc.
unsigned int MYID;
char* PNAME;

/// Blast res
typedef struct blastres {
    uint64_t seqid; 
    uint64_t alignlen; uint64_t qstart; uint64_t qend;
    double evalue; 
    int bitscore;
} BLASTRES;

int BLASTRESSZ = sizeof(BLASTRES);

typedef struct structEvalue {
    BLASTRES* p;
    double evalue;
} STRUCTEVALUE;
string WORKEROUTPUTFILENAME = "";

/// To pass NUMTOTALDBCHUNKS and NUMCOREPERNODE to map() for custom scheduler
typedef struct gf {
    int NUMTOTALDBCHUNKS;
    int NUMCOREPERNODE;
} GF;
GF GF1;
 
/// MR-MPI CALLS
void    set_default_opts(CRef<CBlastOptionsHandle> optsHandle);
void    mr_run_blast(int itask, char *file, KeyValue *kv, void *ptr);
//bool    check_exclusion(string qGi, string sGi, int qCutLocStart, 
                        //uint64_t qCutLocEnd, uint64_t sStart, uint64_t sEnd, 
                        //uint64_t threshold);                    
void    mr_sort_multivalues_by_evalue(char *key, int keybytes, char *multivalue, 
                        int nvalues, int *valuebytes, KeyValue *kv, void *ptr);                         
        
/// TOKENIZER ROUTINES
vector<string>  &split(const string &s, char delim, vector<string> &vecElems);
vector<string>  split(const string &s, char delim);

/// NOTE: NCBI PROVIDES THESE UTILS
string      uint2str(uint64_t number);
uint64_t    str2uint(string str);
int         str2int(string str);

/* -------------------------------------------------------------------------- */
int main(int argc, char **argv)
/* -------------------------------------------------------------------------- */
{
    /// Declare time variables:
    #ifdef DEBUG
    time_t  t0, t1, t2, t3, accMapWTime=0; 
    t0 = time(NULL);
    clock_t c0, c1, c2, c3, accMapCTime=0;     
    c0 = clock();
    #endif
    
    ///
    /// Read "mrblast.conf" to set parameters
    ///
    ConfigFile config("mrblast.conf");
    MAXQUERYNUM = config.read<unsigned int>("MAXQUERYNUM");
    DBCHUNKLISTFILE = config.read<string>("DBCHUNKLISTFILE");
    NUMTOTALDBCHUNKS = config.read<unsigned int>("NUMTOTALDBCHUNKS");          
    GF1.NUMTOTALDBCHUNKS = NUMTOTALDBCHUNKS;
    GF1.NUMCOREPERNODE = config.read<int>("NUMCOREPERNODE");          
    
    /// Read Blast opts
    //string BLASTOPTS = config.read<string>("BLASTOPTS");
    CUTOFFSCORE = config.read<unsigned int>("CUTOFFSCORE");
    EVALUE = config.read<float>("EVALUE");

    ///
    /// MPI setup
    ///
    MPI_Init(&argc, &argv);

    char MPI_procName[MAXSTR];
    int MPI_myId, MPI_nProcs, MPI_length;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myId);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_nProcs);
    MPI_Get_processor_name(MPI_procName, &MPI_length);
    fprintf(stdout, "### INFO: [Rank %d] %s \n", MPI_myId, MPI_procName);
    
    MPI_Barrier(MPI_COMM_WORLD);  
   
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
  
    /// 
    /// Blast search v2
    ///
    MapReduce *mr2 = new MapReduce(MPI_COMM_WORLD);
    
    /*
    * mapstyle = 0 (chunk) or 1 (stride) or 2 (master/slave)
    * all2all = 0 (irregular communication) or 1 (use MPI_Alltoallv)
    * verbosity = 0 (none) or 1 (summary) or 2 (histogrammed)
    * timer = 0 (none) or 1 (summary) or 2 (histogrammed)
    * memsize = N = number of Mbytes per page of memory
    * minpage = N = # of pages to pre-allocate per processor
    * maxpage = N = max # of pages allocatable per processor, 0 = no limit
    * keyalign = N = byte-alignment of keys
    * valuealign = N = byte-alignment of values
    * fpath = string 
    */
    mr2->verbosity = 0;
    mr2->timer = 0;
    mr2->mapstyle = atoi(argv[4]);  /// master/slave mode=2, custom scheduler=3
    mr2->keyalign = sizeof(uint64_t);
    
    MPI_Barrier(MPI_COMM_WORLD);

    uint64_t nvecRes;
    const char *masterFileName = argv[1];
    MYID = MPI_myId;
    PNAME = MPI_procName;
    string prefix(argv[2]);
    ITER = atoi(argv[3]);
    
    ///
    /// Make a file for map() which contains a list of file
    /// names. A file name is a form of "queryFile,dbChunkName"
    /// In map(), the file name is splitted into "queryFile" and
    /// "dbChunkName". We've got 109 DB chunks, FYI.
    ///
    vector<string> vWorkItems;
    vector<string> vDbChunkNames;
    vector<string> vQueryFileNames;
    uint64_t numDbChunks;
    uint64_t numQueryFiles;
    
    if (MPI_myId == 0) {
        const char *dbChunkNameFileName = DBCHUNKLISTFILE.c_str();
        ifstream dbChunkNameFile(dbChunkNameFileName);
        if (!dbChunkNameFile.is_open()) {
            cerr << "### ERROR: dbchunks.txt open error.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        string line;
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
        
        while (!getline(masterFile, line).eof()) {
            vQueryFileNames.push_back(line);
        }        
        masterFile.close();
        
        numQueryFiles = vQueryFileNames.size();
        for (uint64_t i = 0; i < numQueryFiles; ++i) {
            for (uint64_t j = 0; j < numDbChunks; ++j) { 
                vWorkItems.push_back(vQueryFileNames[i] + "," + vDbChunkNames[j]);
            }
        }
        cout << "### INFO: Total number of work items = "
             << vWorkItems.size() << endl;
             
        vQueryFileNames.clear();
        vDbChunkNames.clear();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);  
    
    ///
    /// Now vWorkItems has all the work items in the form of 
    /// (query_file_name,DB_chunk_name)
    ///
    uint64_t numSubMasterFiles = 0;
    uint64_t numWorkItemsForEach = 0;
    
    ///
    /// Split work item master file into sub master files and process the 
    /// work items iteratively while saving the results at each iteration.
    /// eq) if the total number of query files x and ITER = x/2, the master
    /// file is divided into 2 sub master files and processed individually.
    ///
    assert(ITER > 0);
    if (MPI_myId == 0) {
        ///
        /// Divide the whole work items into several work item lists.
        /// Each iteration consists of calling blast and saving output
        /// into an output file.
        ///
        uint64_t numQueryBunch = ITER;
        if (numQueryBunch > numQueryFiles) {
            cerr << "### ERROR: should be numQueryBunch <= numQueryFiles\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        numSubMasterFiles = numQueryFiles / numQueryBunch;
        uint64_t numTotalWorkItems = vWorkItems.size();
        numWorkItemsForEach = numQueryBunch * numDbChunks;
        uint64_t numRemains = numTotalWorkItems % numWorkItemsForEach;
                
        uint64_t k = 0;
        uint64_t i = 0;
        
        for (; i < numSubMasterFiles; ++i) {
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
                numSubMasterFiles++; /// for one remains file
            }
            else {
                cerr << "### ERROR: work item file open error.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        cout << "### INFO: num work item files = " << numSubMasterFiles << endl;
    }
    
    MPI_Barrier(MPI_COMM_WORLD); 

    ///
    /// Iteratively call blast and save results for numSubMasterFiles
    /// times.
    ///

    /// Broadcast the number of sub work item files
    MPI_Bcast(&numSubMasterFiles, 1, MPI_INT, 0, MPI_COMM_WORLD);

    for (uint64_t n = 0; n < numSubMasterFiles; ++n) {
        vector<string> vSubWorkitems;
        string splitMasterFileName = prefix + "-workitems-" + uint2str(n) + ".txt";
        
        ///
        /// This should be gone after remove MapReduce::map(vector<string>...)
        ///
        const char *tempWorkItemFileName = splitMasterFileName.c_str();
        ifstream tempWorkItemFile(tempWorkItemFileName);
        if (!tempWorkItemFile.is_open()) {
            cerr << "### ERROR: split master file open error.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        string line;
        while (!getline(tempWorkItemFile, line).eof()) {
            vSubWorkitems.push_back(line);
        }            
        tempWorkItemFile.close();        
        
        
        ///
        /// map, collate, and reduce
        ///
        #ifdef DEBUG
        c2 = clock();
        t2 = time(NULL);
        #endif
        
        ////////////////////////////////////////////////////////////////////////
        //nvecRes = mr2->map((char*)splitMasterFileName.c_str(), &mr_run_blast, &GF1);
        nvecRes = mr2->map(vSubWorkitems, &mr_run_blast, &GF1);
        ////////////////////////////////////////////////////////////////////////
        
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
        /// New version with sorting done in workders
        /// 1. map(): blast call with assigned work item
        /// 2. collate(): creates KMV with (Qid,DBid) as key
        /// 3. reduce(): convert <(Qid,DBid), {blast_result(s)}> -> <Qid, {blast_result(s)}>
        /// 4. collate(): creates KMV with Qid as key
        /// 5. reduce(): sort multivalues by bit score
        /// 6. gather(): gather KVs to master
        /// 7. print(): printing
        ///
        ////////////////////////////////////////////////////////////////////////
        mr2->collate(NULL);
        WORKEROUTPUTFILENAME = prefix + "-" + uint2str(n) + "-" 
                               + uint2str(MPI_myId) + ".txt";
        uint64_t nunique = mr2->reduce(&mr_sort_multivalues_by_evalue, NULL);
        ////////////////////////////////////////////////////////////////////////
        
        ///
        /// Save history
        ///
        if (MPI_myId == 0) {                
            string histFileName = prefix + "-history.txt";
            FILE *histFile = fopen(histFileName.c_str(), "a");
            time_t timer;
            timer = time(NULL);
            fprintf(histFile, "%s,%s", splitMasterFileName.c_str(), asctime(localtime(&timer)));
            fclose(histFile);
        }
            
        MPI_Barrier(MPI_COMM_WORLD);           
    }
     
    delete mr2;
    delete pTargetDb;   
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

void mr_run_blast(int itask, char *file, KeyValue *kv, void *ptr)
{
    #ifdef DEBUG
    time_t t0, t1, t2, t3, accBlastCallTime=0;
    time_t t4, t5, accDBLoadingTime=0;
    t0 = time(NULL);
    #endif 
    
    /// 
    /// Make a option handle and cblastinputsource
    ///
    EProgram program = ProgramNameToEnum("blastn");
    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));
    set_default_opts(opts);
    opts->Validate();

    /// DEBUG
    if (MYID == 0)
        opts->GetOptions().DebugDumpText(cerr, "opts", 1);
    ///

    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
        throw std::runtime_error("Could not initialize object manager");
    }

    bool isProtein = false;
    SDataLoaderConfig dlconfig(isProtein);
    CBlastInputSourceConfig iconfig(dlconfig);
    
    //string qFileName = "test2.query";
    //ifstream queryFile(qFileName.c_str());
    //CBlastFastaInputSource fasta_input(qFileName, iconfig);

    ///
    /// Split work item => query file + db chunk name
    ///
    string workItem(file);
    vector<string> vWorkItem = split(workItem, ',');
    ifstream queryFile(vWorkItem[0].c_str());
    if (!queryFile.is_open()) {
        cerr << "### ERROR: queryFile open error.\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

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
    if(pTargetDb == 0 || dbChunkName != prevDbChunkName) {
        delete pTargetDb;
        pTargetDb = new CSearchDatabase(dbChunkName,
                                        CSearchDatabase::eBlastDbIsNucleotide);
    }
    prevDbChunkName = dbChunkName;
    
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
                    Seq read error.\n", MYID, PNAME, itask, file);
            MPI_Abort(MPI_COMM_WORLD, 1);       
            exit(0);
        }
        
        /// NOTE: This should be improved!
        string query = header + '\n' + seq;
        vHeaders.push_back(header);
        ///

        ///
        /// ACCUMMULATE SEQS WHEN MAXQUERYNUM > 1
        ///
        uint64_t j = 1;
        while (j < MAXQUERYNUM && !getline(queryFile, header).eof()) {
            getline(queryFile, seq);
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
        CLocalBlast blaster(queryFactory, opts, *pTargetDb);
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

                    /// qNum is not unique. It's internal qid in fasta_input
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

                    uint64_t alignLen, qStart, qEnd;
                    //uint64_t sStart, sEnd;
                    //CRange<TSeqPos> qRange, sRange;
                    //qRange    = s.GetSeqRange(QUERY);
                    //sRange    = s.GetSeqRange(SUBJECT);
                    qStart      = s.GetSeqStart(QUERY);
                    qEnd        = s.GetSeqStop(QUERY);
                    //sStart      = s.GetSeqStart(SUBJECT);
                    //sEnd        = s.GetSeqStop(SUBJECT);
                    alignLen    = s.GetAlignLength();

                    double eValue;
                    //int genericScore; 
                    int bitScore;
                    s.GetNamedScore(CSeq_align::eScore_EValue, eValue);
                    //s.GetNamedScore(CSeq_align::eScore_Score, genericScore);
                    s.GetNamedScore(CSeq_align::eScore_BitScore, bitScore);
                    
                    ///
                    /// TOKENIZE QUERY HEADER
                    /// 
                    string qHeader = vHeaders[str2uint(qNum)-1];
                    vector<string> vQeuryId  = split(qHeader, '|');
                    //string qGi               = vQeuryId[1];             /// GI
                    string qid               = vQeuryId[2];             /// QUERY ID
                    ////uint64_t origLen       = str2uint(vQeuryId[3]);   /// LENGTH OF THE ORIG SEQ
                    //int qCutLocStart         = str2int(vQeuryId[4]);    /// CUT COORDINATES - START
                    //uint64_t qCutLocEnd      = str2uint(vQeuryId[5]);   /// CUT COORDINATES - END

                    /// 
                    /// ADD A CSV BLAST RESULT TO KV
                    /// 
                    //if (!check_exclusion(qGi, seqId, qCutLocStart, qCutLocEnd, sStart, sEnd, EXCLUSION)) {
                        //char blastRes[MAXSTR];

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
                        //sprintf(blastRes, "%s,%s,%ld,%ld,%ld,%1.e,%d",
                                //qHeader.substr(0, MAXSTR).c_str(), seqId.c_str(),
                                //alignLen, qStart, qEnd,
                                //eValue, bitScore);
                        
                        ///
                        /// 12.20.2010
                        /// 
                        BLASTRES r;
                        //r.qheader = qHeader.substr(0, MAXSTR);
                        r.seqid = str2uint(seqId);
                        r.alignlen = alignLen;
                        r.qstart = qStart;
                        r.qend = qEnd;
                        r.evalue = eValue;
                        r.bitscore = bitScore;     
                        //char *newKey = (char*)((qid).c_str()); 
                        uint64_t newKey = str2uint(qid);
                        
                        /// 
                        /// ADD <KEY = "QUERYID,CHUNKNAME", VALUE="BLASTvecResULT">
                        /// TO KV
                        /// 
                        //kv->add(newKey, strlen(newKey) + 1, blastRes, strlen(blastRes) + 1);
                        //kv->add((char*)&newKey, sizeof(uint64_t), blastRes, strlen(blastRes) + 1);
                        kv->add((char*)&newKey, sizeof(uint64_t), (char*)&r, BLASTRESSZ);
                        
                    //}
                    //else {
                        //fprintf(stdout, "### INFO: EXCLUSION - qGi=%s sGi=%s \
                                //qCutLocStart=%ld, qCutLocEnd=%ld, sStart=%ld, sEnd=%ld\n",
                                //qGi.c_str(), seqId.c_str(), qCutLocStart, qCutLocEnd, 
                                //sStart, sEnd);

                        ////if (exclusionSaveFile.is_open()) {
                        ////char excluded[MAXSTR];
                        ////sprintf(stderr, "qGi=%s sGi=%s qCutLocStart=%d, qCutLocEnd=%d, sStart=%d, sEnd=%d\n", qGi.c_str(), seqId.c_str(), qCutLocStart, qCutLocEnd, sStart, sEnd);
                    //}
                }
            }
        }
    } /// END WHILE
    
    #ifdef DEBUG
    t1 = time(NULL);
    /// FORMAT: rank,total_runblast2,blaster.Run,DBLoadingTime
    printf ("### TIME (in runblast2): rank,total_runblast2,blaster.Run,DBLoadingTime,%d,%ld,%ld,%ld\n", 
            MYID, (long)(t1 - t0), (long)accBlastCallTime, (long)accDBLoadingTime);
    #endif
}


/** Sort function - Passed to MR-MPI sort_values() for sorting blast result
 * string by bit score.
 * @param e1
 * @param e2
 */
 
bool mycompare(STRUCTEVALUE e1, STRUCTEVALUE e2)
{
    return (e1.evalue < e2.evalue);
}

/** Sort by wvalue - Passed to MR-MPI reduce() for sorting KMVs by evalue.
 * @param key
 * @param keybytes
 * @param multivalue - collected blast result strings. There could be more than
 * nvalues results in it. The separator '>' should be used for splitting the 
 * resutls.
 * @param nvalues
 * @param valuebytes
 * @param kv
 * @param ptr
 */
 
void mr_sort_multivalues_by_evalue(char *key, int keybytes, char *multivalue,
                                  int nvalues, int *valuebytes, KeyValue *kv, 
                                  void *ptr) 
{
    /// Check if there is KMV overflow
    assert(multivalue != NULL && nvalues != 0);
    
    ///
    /// Make STRUCTEVALUE = {BLASTRES* p; double evalue;}
    /// and sort
    ///
    vector<STRUCTEVALUE> vforsort;
    for (size_t n = 0; n < nvalues; n++) {                
        BLASTRES* res = (BLASTRES*)multivalue;
        STRUCTEVALUE structEvalue;
        structEvalue.p = res;
        structEvalue.evalue = res->evalue;
        vforsort.push_back(structEvalue);
        multivalue += BLASTRESSZ;         
    }
    sort(vforsort.begin(), vforsort.end(), mycompare);
    
    ///
    /// After sorting, each worker saves the set of results into a file.
    ///
    FILE* fp = fopen((char*)WORKEROUTPUTFILENAME.c_str(), "a");
    for (size_t n = 0; n < nvalues; n++) {
        BLASTRES* res = (BLASTRES*)(vforsort[n].p);
        //cout << *((uint64_t *)key) << "," << res->seqid << "," << res->evalue << endl;
        fprintf(fp, "%ld,%ld,%ld,%ld,%ld,%1.e,%d\n", *((uint64_t *)key),
                res->seqid, res->alignlen, res->qstart, res->qend,
                res->evalue, res->bitscore);
    }
    cout << "### INFO: file saved, " << WORKEROUTPUTFILENAME << endl; 
    fclose(fp);
    vforsort.clear();
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
//bool check_exclusion(string qGi, string sGi, int qCutLocStart,
                     //uint64_t qCutLocEnd, uint64_t sStart, uint64_t sEnd,
                     //uint64_t threshold)
//{
    ///// 
    ///// To exclude Blast result from the original sequence from which
    ///// the input query is originated (sampled). Basically if qGi == sGi
    ///// and qCutLocStart is similar with sStart and qCutLocEnd is similar
    ///// with sEnd in terms of coordinates, the result should be excluded.
    ///// 
    ///// Orig seq: ----------------XXXXXXXXXXXXXXXX----------------------
    /////                           |              |
    /////                     qCutLocStart      qCutLocEnd
    ///// 
    ///// Query:                    XXXXXXXXXXXXXXXX
    /////                              |          |
    /////                           qStart       qEnd
    ///// 
    ///// Subject:   ------------------XXXXXXXXXXXX-----------------------
    /////                              |          |
    /////                           sStart      sEnd
    ///// 
    //bool ret = false;

    //if (qGi == sGi) {
        //if (qCutLocStart < 0) {
            ///// 
            ///// In >gi|222299657|18|3605|-400|3604
            ///// -400|3604 means query[-400:3604] in Python.
            ///// 
            //qCutLocStart = qCutLocEnd + 1 - qCutLocStart;
            //qCutLocEnd += 1;
        //}
        //if ((qCutLocStart - threshold <= sStart 
            //&& sStart <= qCutLocStart + threshold) 
            //&& (qCutLocEnd - threshold <= sEnd 
            //&& sEnd <= qCutLocEnd + threshold))
            //ret = true;
    //}

    //return ret;
//}

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

 

/// EOF

