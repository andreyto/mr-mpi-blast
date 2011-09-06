//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

////////////////////////////////////////////////////////////////////////////////
//
//  MR-MPI-BLAST: Parallelizing BLAST using MapReduce-MPI
//
//  Author: Seung-Jin Sul
//         (ssul@jcvi.org)
//
//  Last updated: 08/01/2011
//
////////////////////////////////////////////////////////////////////////////////

#include "mrblast.hpp"


void mrmpi_blast(); 
void mr_map_run_blast(int itask, KeyValue *kv, void *ptr);
inline void mpi_collect_node_name(int rank, int nProcs, MPI_Comm mpiComm);
inline void mr_reduce_sort_multivalues_by_evalue(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
inline void mr_reduce_save_hits(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
inline bool compare_evalue(structEValue_t e1, structEValue_t e2);


class CMrMpiBlastApplication : public CNcbiApplication
{
private:
    virtual void Init(void);
    virtual int  Run(void); 
};

void CMrMpiBlastApplication::Init(void)
{
    ///
    /// Read conf file, mrblast.ini and set parameters
    ///
    ifstream config(CONF_FILE_NAME.c_str(), ios::in);
    if (!config) {
        MPI_Abort(MPI_COMM_WORLD, 1);
        cout << "ERROR: failed to open configuration file (mrblast.ini)" << endl;
        exit(1);
    }

    set<string> options;
    mapSS_t parameters;
    options.insert("*");

    try {
        for (pod::config_file_iterator i(config, options), e ; i != e; ++i) {
            parameters[i->string_key] = i->value[0];
        }
        try {
            g_verbosity      = boost::lexical_cast<int>(parameters["VERBOSITY"]);
            g_timer          = boost::lexical_cast<int>(parameters["TIMER"]);
            g_memSize        = boost::lexical_cast<int>(parameters["MEMSIZE"]);
            g_outOfCore      = boost::lexical_cast<int>(parameters["OUTOFCORE"]);

            g_bLogEnabled    = boost::lexical_cast<bool>(parameters["LOGENABLED"]);
            g_timingEnabled  = boost::lexical_cast<int>(parameters["TIMING"]);
            g_optDumpEnabled = boost::lexical_cast<int>(parameters["OPTDUMP"]);
                        
            g_identCutoff    = boost::lexical_cast<float>(parameters["IDENTCUTOFF"]);
            g_coverCutoff    = boost::lexical_cast<float>(parameters["COVERCUTOFF"]);
            
            g_logFileName    = parameters["LOGFNAME"];
            g_queryFileName  = parameters["QUERYFILENAME"];
            g_indexFileName  = parameters["INDEXFILENAME"];
            g_dbFileName     = parameters["DBLISTFILENAME"];
            g_outFilePrefix  = parameters["OUTFILEPREFIX"];
            
            g_blockSize      = boost::lexical_cast<uint32_t>(parameters["BLOCKSIZE"]);
            g_mapStyle       = boost::lexical_cast<int>(parameters["MAPSTYLE"]);
            g_nIter          = boost::lexical_cast<int>(parameters["NUMITER"]);
            g_bIsProtein     = boost::lexical_cast<bool>(parameters["ISPROTEIN"]);
        }
        catch (const boost::bad_lexical_cast &) {
            cout << "Exception: bad_lexical_cast" << endl;
            cout.flush();
        }
    }
    catch (exception& e) {
        cout << "Exception: " << e.what() << endl;
        cout.flush();
    }
     
    if (g_bIsProtein) g_cmdLineArgs.Reset(new CBlastpAppArgs());
    else g_cmdLineArgs.Reset(new CBlastnAppArgs());
    HideStdArgs(fHideLogfile | fHideConffile | fHideFullVersion | 
                fHideXmlHelp | fHideDryRun);
    SetupArgDescriptions(g_cmdLineArgs->SetCommandLine());  
}

int CMrMpiBlastApplication::Run(void) 
{     
    /// Allow the fasta reader to complain on invalid sequence input
    SetDiagPostLevel(eDiag_Warning);
    
    const CArgs& args = GetArgs();
    string allArgs;
    RecoverSearchStrategy(args, g_cmdLineArgs);
    CRef<CBlastOptionsHandle> optsHndl(&*g_cmdLineArgs->SetOptions(args));
    g_optsHndl = optsHndl;
    
    ///
    /// Collect MPI node names and rank number for mapstyle=3 scheduler
    ///
    mpi_collect_node_name(g_MPI_worldRank, g_MPI_nProcs, MPI_COMM_WORLD);
    
    if (g_MPI_worldRank == 0) {
        for (multimapSI_t::const_iterator iter = g_multimapProcNameRank.begin();
             iter != g_multimapProcNameRank.end(); ++iter )
            cout << iter->first << '\t' << iter->second << '\n';
    }
   
    ///
    /// Creat memory-mapped file for query
    ///
    g_realFileSize = boost::filesystem::file_size(g_queryFileName);
    g_memmapQueryFile.open(g_queryFileName, g_realFileSize, 0);
    if (!g_memmapQueryFile.is_open()) {
        MPI_Abort(MPI_COMM_WORLD, 1);
        cout << "ERROR: failed to create mmap query file\n";
        exit(1);
    }
 
    ///
    /// Load DB list 
    ///
    string line;
    ifstream dbListFile(g_dbFileName.c_str(), ios::in);
    if (!dbListFile.is_open()) {        
        MPI_Abort(MPI_COMM_WORLD, 1);
        cout << "ERROR: failed to open DB list file.\n";
        exit(1);
    }
    while (!getline(dbListFile, line).eof() && line.length() > 0)
        g_vecDbFile.push_back(line);
    dbListFile.close();
    g_nDbFiles = g_vecDbFile.size();

    ///
    /// Load index file
    ///
    ifstream indexFile(g_indexFileName.c_str(), ios::in);
    if (!indexFile.is_open()) {
        MPI_Abort(MPI_COMM_WORLD, 1);
        cout << "ERROR: failed to open index file.\n";
        exit(1);
    }
    while (!getline(indexFile, line).eof() && line.length() > 0) {
        vector<string> tok;
        boost::split(tok, line, boost::is_any_of(","));
        assert(tok.size() == 2);
        
        structQueryIndex_t queryIndex;
        queryIndex.qStart  = boost::lexical_cast<uintmax_t>(tok[0]);
        queryIndex.qLength = boost::lexical_cast<uint32_t>(tok[1]);
        g_vecQueryIndex.push_back(queryIndex);
    }
    indexFile.close();
    g_nQueries = g_vecQueryIndex.size();

    ///
    /// Fill g_vQueryStartLoc vector from the index file using g_blockSize
    ///
    uint32_t qSizeCurr = g_blockSize;
    for (size_t q = 0; q < g_nQueries; ++q) {
        if (qSizeCurr >= g_blockSize) {
            g_vecBlockBeginLoc.push_back(g_vecQueryIndex[q].qStart);
            qSizeCurr = 0;
        }
        qSizeCurr += g_vecQueryIndex[q].qLength;
    }
    g_nQueryBlocks = g_vecBlockBeginLoc.size();

    ///
    /// Create work items
    /// Here we iterate query block ID with fixing each DB partition name
    ///
    for (size_t d = 0; d < g_nDbFiles; ++d) {
        for (size_t b = 0; b < g_nQueryBlocks - 1; ++b) {
            structWorkItem_t aWorkItem;
            aWorkItem.blockBegin = g_vecBlockBeginLoc[b];
            aWorkItem.blockEnd   = g_vecBlockBeginLoc[b + 1] - 1;
            aWorkItem.dbNo       = d;
            g_vecWorkItem.push_back(aWorkItem);
        }
        
        structWorkItem_t aWorkItem;
        aWorkItem.blockBegin = g_vecBlockBeginLoc[g_nQueryBlocks - 1];
        aWorkItem.blockEnd   = g_realFileSize;
        aWorkItem.dbNo       = d;
        g_vecWorkItem.push_back(aWorkItem);
    }
    g_nWorkItems = g_vecWorkItem.size();
    
    ///
    /// Save the list of workitem just in case
    ///
    if (g_MPI_worldRank == 0) {
        string workItemFName = g_outFilePrefix + "-workitems.txt";
        ofstream workItemFile(workItemFName.c_str(), ios::out);
        
        for (size_t w = 0; w < g_nWorkItems; ++w)
            workItemFile << g_vecWorkItem[w].blockBegin << ","
                         << g_vecWorkItem[w].blockEnd << ","
                         << g_vecWorkItem[w].dbNo << endl;
        workItemFile.close();
    }        
    
    if (g_MPI_worldRank == 0) {
        cout << "Input query file = " << g_queryFileName << endl;
        cout << "Input query index file = " << g_indexFileName << endl;
        cout << "Database list file = " << g_dbFileName << endl;
        cout << "Block size (bp) = " << g_blockSize << endl;
        cout << "Number of query blocks = " << g_nQueryBlocks << endl;
        cout << "Number of DB files = " << g_nDbFiles << endl;
        cout << "Number of total work items = " << g_nWorkItems << endl;
        cout << "Number of iterations = " << g_nIter << endl;
        cout << "Identity cutoff = " << g_identCutoff << endl;
        cout << "Coverage cutoff = " << g_coverCutoff << endl;        
        cout.flush();
    }

    ///
    /// Calculate sub work item size for multiple iterations
    ///
    if (g_nIter != 1) {
        nSubWorkItemSets = g_nIter;
        nWorkItemsPerIter = (g_nQueryBlocks / g_nIter) * g_nDbFiles;            
        nRemains = (g_nQueryBlocks % g_nIter) * g_nDbFiles;
    }
    else { 
        nWorkItemsPerIter = g_nWorkItems;
        nRemains = 0;
        nSubWorkItemSets = 1;
    }
    /// for one more added file of remaining work items
    if (nRemains) ++nSubWorkItemSets; 

    if (g_MPI_worldRank == 0 && g_nIter != 1) {        
        cout << "Number of sub work item sets = " << nSubWorkItemSets << endl;
        cout << "Number of work items per iteration  = " << nWorkItemsPerIter << endl;
        cout << "Number of work items remaining = " << nRemains << endl;
        cout.flush();
    }       
        
    MPI_Barrier(MPI_COMM_WORLD);

    ///
    /// Run mr-mpi calls
    ///
    mrmpi_blast(); 
    if (g_MPI_worldRank == 0) cout << "Done!" << endl;
    
    ///
    /// Clean up
    ///
    g_vecDbFile.clear();
    g_vecQueryIndex.clear();
    g_vecBlockBeginLoc.clear();
    g_vecWorkItem.clear();
    
    if (!g_searchDatabase.IsNull()) g_searchDatabase.Release();
    if (g_bLogEnabled || g_timingEnabled) g_logFileStream.close();
    if (g_memmapQueryFile.is_open()) g_memmapQueryFile.close(); /// close mmapped file  
    

    return 0;
}
 
 
int main(int argc, char** argv)
{
    ///
    /// MPI setup
    ///
    int MPI_procNameLen;    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_MPI_worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &g_MPI_nProcs);
    MPI_Get_processor_name(g_MPI_procName, &MPI_procNameLen);
        
    /// Execute main application function
    int ret = CMrMpiBlastApplication().AppMain(argc, argv);    
    
    MPI_Finalize();   
    
    return ret;
}


/** MapReduce fuction for BLAST search
 */

void mrmpi_blast()
{
    int rank = g_MPI_worldRank;
    
    ///
    /// Log file init
    ///
    if (g_bLogEnabled || g_timingEnabled) {
        g_logFileName = g_outFilePrefix + "-" + boost::lexical_cast<string>(rank) 
                        + "-" + g_logFileName;
        g_logFileStream.open(g_logFileName.c_str(), ios::out);
        g_pLogStream = &g_logFileStream;
    }
    g_logMsg = "Rank:" + boost::lexical_cast<string>(rank) + " ";

    double profileTime;
    struct timeval wallTimeStart;
    struct timeval wallTimeEnd;
    struct timeval userTimeStart;
    struct timeval userTimeEnd;
    struct timeval sysTimeStart;
    struct timeval sysTimeEnd;
    struct rusage ruTotal;

    if (g_timingEnabled) {
        profileTime = MPI_Wtime(); /// MPI_Wtime
        gettimeofday(&wallTimeStart, NULL); /// Wall-clock time
        getrusage(RUSAGE_SELF, &ruTotal); /// Process time
        userTimeStart = ruTotal.ru_utime;
        sysTimeStart = ruTotal.ru_stime;
        LOG << "mr-mpi-blast starts,"
            << profileTime << ","
            << wallTimeStart.tv_sec * 1000000 + wallTimeStart.tv_usec << ","
            << userTimeStart.tv_sec * 1000000 + userTimeStart.tv_usec << ","
            << sysTimeStart.tv_sec * 1000000 + sysTimeStart.tv_usec
            << endl;
        LOG.flush();
    }

    ///
    ///
    /// MR-MPI init (MAY2011 version)
    ///
    /// mapstyle    = 0 (chunk) or 1 (stride) or 2 (master/slave)
    /// all2all     = 0 (irregular communication) or 1 (use MPI_Alltoallv)
    /// verbosity   = 0 (none) or 1 (summary) or 2 (histogrammed)
    /// timer       = 0 (none) or 1 (summary) or 2 (histogrammed)
    /// memsize     = N = number of Mbytes per page of memory
    /// minpage     = N = # of pages to pre-allocate per processor
    /// maxpage     = N = max # of pages allocatable per processor
    /// freepage    = 1 if memory pages are freed in between operations, 0 if held
    /// outofcore   = 1 if even 1-page data sets are forced to disk, 0 if not, -1
    ///               if cannot write to disk
    /// zeropage    = 1 if zero out every allocated page, 0 if not
    /// keyalign    = N = byte-alignment of keys
    /// valuealign  = N = byte-alignment of values
    /// fpath       = string
    ///
    MapReduce *pMr = new MapReduce(MPI_COMM_WORLD);
    pMr->verbosity = g_verbosity;
    pMr->timer     = g_timer;
    pMr->memsize   = g_memSize;
    pMr->keyalign  = sizeof(uint32_t);  /// The key is query GI
    pMr->mapstyle  = g_mapStyle;        /// master/slave=2, custom scheduler=3
    pMr->outofcore = g_outOfCore;       /// "-1" to disable out-of-core operation
    MPI_Barrier(MPI_COMM_WORLD);

    for (size_t iter = 0; iter < nSubWorkItemSets; ++iter) {
        if (g_MPI_worldRank == 0) cout << "Iteration: " << iter << endl;
        if (g_bLogEnabled) { 
            LOG << g_logMsg << "Iteration: " << iter << endl;
            LOG.flush();
        }
        
        g_currIter = iter;
        
        ///
        /// map, collate, reduce
        ///
        double mapTime;
        if (g_bLogEnabled) {
            mapTime = MPI_Wtime();
            LOG << g_logMsg << "map() starts: " <<  mapTime << endl;
            LOG.flush();
        }
 
        /////////////////////////////////////////////////////////////
        if (g_nIter != 1)
            if (nRemains != 0 && iter == nSubWorkItemSets - 1)
                pMr->map(nRemains, &mr_map_run_blast, NULL);
            else 
                pMr->map(nWorkItemsPerIter, &mr_map_run_blast, NULL);
        else 
            pMr->map(g_nWorkItems, &mr_map_run_blast, NULL);
        /////////////////////////////////////////////////////////////
                
        if (g_bLogEnabled) {
            LOG << g_logMsg << "map() ends: " <<  MPI_Wtime() - mapTime << endl;
            LOG.flush();
        }

        double collateTime;
        if (g_bLogEnabled) {
            collateTime = MPI_Wtime();
            LOG << g_logMsg << "collate() starts: " << collateTime << endl;
            LOG.flush();
        }

        ///////////////////
        pMr->collate(NULL);
        ///////////////////

        if (g_bLogEnabled) {
            LOG << g_logMsg << "collate() ends: " << MPI_Wtime() - collateTime 
                << endl;
            LOG.flush();
        }

        double reduceTime;
        if (g_bLogEnabled) {
            reduceTime = MPI_Wtime();
            LOG << g_logMsg << "reduce() starts: " << reduceTime << endl;
            LOG.flush();
        }
        
        g_hitFileName = g_outFilePrefix + "-hits-" 
                    + boost::lexical_cast<string>(iter) + "-"
                    + boost::lexical_cast<string>(rank) + ".txt";
                        
        ///////////////////////////////////////////////////
        //pMr->reduce(&mr_reduce_sort_multivalues_by_evalue, NULL);
        pMr->reduce(&mr_reduce_save_hits, NULL);
        ///////////////////////////////////////////////////

        if (g_bLogEnabled) {
            LOG << g_logMsg << "reduce() ends: " <<  MPI_Wtime() - reduceTime 
                << endl;
            LOG.flush();
        }
       
        if (g_timingEnabled) {
            /// Wall-clock time
            gettimeofday(&wallTimeEnd, NULL);
            double tS = wallTimeStart.tv_sec * 1000000 + (wallTimeStart.tv_usec);
            double tE = wallTimeEnd.tv_sec * 1000000  + (wallTimeEnd.tv_usec);

            /// process time
            getrusage(RUSAGE_SELF, &ruTotal);
            userTimeEnd = ruTotal.ru_utime;
            sysTimeEnd = ruTotal.ru_stime;
            double tS_user = userTimeStart.tv_sec * 1000000
                             + (userTimeStart.tv_usec);
            double tE_user = userTimeEnd.tv_sec * 1000000
                             + (userTimeEnd.tv_usec);
            double tS_sys  = sysTimeStart.tv_sec * 1000000
                             + (sysTimeStart.tv_usec);
            double tE_sys  = sysTimeEnd.tv_sec * 1000000
                             + (sysTimeEnd.tv_usec);

            double profile_etime = MPI_Wtime();
            LOG << "mr-mpi-blast ends,"
                << profile_etime << ","
                << wallTimeEnd.tv_sec * 1000000 + wallTimeEnd.tv_usec << ","
                << userTimeEnd.tv_sec * 1000000 + userTimeEnd.tv_usec << ","
                << sysTimeEnd.tv_sec * 1000000 + sysTimeEnd.tv_usec << endl;

            /// MP_Wtime
            profileTime = profile_etime - profileTime;
            LOG << "Total MPI_Wtime time," << profileTime << endl;
            LOG << "Total wall-clock time,"
                << tE - tS << "," << (tE - tS) / 1000000 << endl;
            LOG << "Total process time (user),"
                << tE_user - tS_user << "," << (tE_user - tS_user) / 1000000 
                << endl;
            LOG << "Total process time (sys),"
                << tE_sys - tS_sys << "," << (tE_sys - tS_sys) / 1000000 << endl;
            LOG << "Total total process time (user+sys),"
                << (tE_user - tS_user) + (tE_sys - tS_sys) << ","
                << ((tE_user - tS_user) + (tE_sys - tS_sys)) / 1000000 << endl;
            LOG.flush();
        }
    }
    
    delete pMr; /// MR-MPI
}


/** MR-MPI Map function - settup NCBI C++ Toolkit env and call blast
 * @param itask
 * @param kv
 * @param ptr
 */
 
void mr_map_run_blast(int itask,
                  KeyValue *kv,
                  void *ptr)
{
    int rank = g_MPI_worldRank;
    int iter = g_currIter;
           
    uint32_t dbno;
    uintmax_t qBlockStart, qBlockEnd;    
    if (g_nIter > 1) {
        dbno = g_vecWorkItem[itask + nWorkItemsPerIter * iter].dbNo;
        qBlockStart = g_vecWorkItem[itask + nWorkItemsPerIter * iter].blockBegin;
        qBlockEnd = g_vecWorkItem[itask + nWorkItemsPerIter * iter].blockEnd;
    }
    else {
        dbno = g_vecWorkItem[itask].dbNo;
        qBlockStart = g_vecWorkItem[itask].blockBegin;
        qBlockEnd   = g_vecWorkItem[itask].blockEnd;
    }

    struct timeval blastcallStartTime;
    struct timeval blastcallEndTime;
    struct timeval blastcallStart_u_Time;
    struct timeval blastcallEnd_u_Time;
    struct timeval blastcallStart_s_Time;
    struct timeval blastcallEnd_s_Time;
    struct rusage  ru_blastcall;
    struct timeval qBuildStartTime;
    struct timeval qBuildStart_u_Time;
    struct timeval qBuildStart_s_Time;
    struct rusage  ru_qBuild;    
    struct timeval dbLoadingStartTime;
    struct timeval dbLoadingStart_u_Time;
    struct timeval dbLoadingStart_s_Time;
    struct rusage  ru_dbLoading;
     
    if (rank == 1 && g_optDumpEnabled == 1) {
        g_optDumpEnabled = 0;
        string strategyFName = g_outFilePrefix + "-search_strategy.txt";
        ofstream strategyOutFile(strategyFName.c_str(), ios::out);
        g_optsHndl->GetOptions().DebugDumpText(strategyOutFile, "optsHndl", 1);
        strategyOutFile.close();
    }
    
    string dbFileName = g_vecDbFile[dbno];
     
    ///
    /// Read a block of sequeces from qBlockStart
    ///
    double query_build_time;
    if (g_timingEnabled) {
        ++g_mapCallNo;
        query_build_time = MPI_Wtime();
        gettimeofday(&qBuildStartTime, NULL); /// Wall-clock time
        getrusage(RUSAGE_SELF, &ru_qBuild);   /// Process time
        qBuildStart_u_Time = ru_qBuild.ru_utime;
        qBuildStart_s_Time = ru_qBuild.ru_stime;
        LOG << g_logMsg << "query build starts,"
            << query_build_time << ","
            << qBuildStartTime.tv_sec * 1000000
             + qBuildStartTime.tv_usec << ","
            << qBuildStart_u_Time.tv_sec * 1000000
             + qBuildStart_u_Time.tv_usec << ","
            << qBuildStart_s_Time.tv_sec * 1000000
             + qBuildStart_s_Time.tv_usec << ","
            << dbFileName << "," << g_mapCallNo << "," << g_MPI_procName
            << "," << qBlockStart << endl;
        LOG.flush();
    }

    assert(g_memmapQueryFile.is_open());
    const char *pMmapQueryFile = (char*)g_memmapQueryFile.data();
    ////////////////////////////////////////////////////////////////////
    string query(pMmapQueryFile + qBlockStart, qBlockEnd - qBlockStart);    
    ////////////////////////////////////////////////////////////////////
    assert(query.length() > 0);
    assert(query[0] == '>');

    int status = BLAST_EXIT_SUCCESS;
    try {
        const CBlastOptions& opt = g_optsHndl->GetOptions();
        
        /*** Get the query sequence(s) ***/
        CRef<CQueryOptionsArgs> query_opts = g_cmdLineArgs->GetQueryOptionsArgs();
        SDataLoaderConfig dlconfig(query_opts->QueryIsProtein());
        dlconfig.OptimizeForWholeLargeSequenceRetrieval();
        CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
                                     query_opts->UseLowercaseMasks(),
                                     query_opts->GetParseDeflines(),
                                     query_opts->GetRange(),
                                     !g_cmdLineArgs->ExecuteRemotely());
        iconfig.SetLowercaseMask(true); /// Enforce lowercase mask option again
        CBlastFastaInputSource fasta(query, iconfig);
        CBlastInput input(&fasta);

        ///
        /// Target db name setting
        ///    
        double db_loading_time;
        if (g_timingEnabled) {
            ++g_mapCallNo;
            db_loading_time = MPI_Wtime();
            gettimeofday(&dbLoadingStartTime, NULL); /// Wall-clock time
            getrusage(RUSAGE_SELF, &ru_dbLoading);   /// Process time
            dbLoadingStart_u_Time = ru_dbLoading.ru_utime;
            dbLoadingStart_s_Time = ru_dbLoading.ru_stime;
            LOG << g_logMsg << "db_loading starts,"
                << db_loading_time << ","
                << dbLoadingStartTime.tv_sec * 1000000
                 + dbLoadingStartTime.tv_usec << ","
                << dbLoadingStart_u_Time.tv_sec * 1000000
                 + dbLoadingStart_u_Time.tv_usec << ","
                << dbLoadingStart_s_Time.tv_sec * 1000000
                 + dbLoadingStart_s_Time.tv_usec << ","
                << dbFileName << "," << g_mapCallNo << "," << g_MPI_procName
                << "," << qBlockStart << endl;
            LOG.flush();
        }
        
        if (g_searchDatabase.IsNull() || dbFileName != g_prevDbName) {
            if (!g_searchDatabase.IsNull()) g_searchDatabase.Release();
            if (g_bIsProtein) {
                CRef<CSearchDatabase> searchDatabase(new CSearchDatabase(
                    dbFileName, CSearchDatabase::eBlastDbIsProtein));
                g_searchDatabase = searchDatabase;
            }
            else {
                CRef<CSearchDatabase> searchDatabase(new CSearchDatabase(
                    dbFileName, CSearchDatabase::eBlastDbIsNucleotide));
                g_searchDatabase = searchDatabase;
            }
        }
        g_prevDbName = dbFileName;
            
        /*** Initialize the database/subject ***/
        CRef<CBlastDatabaseArgs> db_args(g_cmdLineArgs->GetBlastDatabaseArgs());
        CRef<CLocalDbAdapter> db_adapter;
        CRef<CScope> scope;
        db_adapter.Reset();
        CRef<CSearchDatabase> search_db = g_searchDatabase;
        if (scope.Empty()) {
            scope.Reset(new CScope(*CObjectManager::GetInstance()));
        }
        _ASSERT(scope.NotEmpty());
        _ASSERT(search_db.NotEmpty());
        
        try { 
            /// Try to open the BLAST database even for remote searches, as if
            /// it is available locally, it will be better to fetch the
            /// sequence data for formatting from this (local) source
            CRef<CSeqDB> seqdb = search_db->GetSeqDb();
            db_adapter.Reset(new CLocalDbAdapter(*search_db));
            scope->AddDataLoader(RegisterOMDataLoader(seqdb));
        } catch (const CSeqDBException&) {
                /// The BLAST database couldn't be found, report this for local
                /// searches.
        }
        _ASSERT(db_adapter && scope);
            
        if (opt.GetUseIndex() && !g_cmdLineArgs->ExecuteRemotely()) {
            BlastSeqSrc* seqsrc = db_adapter->MakeSeqSrc();
            CRef<CBlastOptions> my_options(&(g_optsHndl->SetOptions()));
            CSetupFactory::InitializeMegablastDbIndex(seqsrc, my_options);
        }
        
        /*** Get the formatting options ***/
        CRef<CFormattingArgs> fmt_args(g_cmdLineArgs->GetFormattingArgs());
        CBlastFormat formatter(opt, 
                               *db_adapter,
                               fmt_args->GetFormattedOutputChoice(),
                               query_opts->GetParseDeflines(),
                               g_cmdLineArgs->GetOutputStream(),
                               fmt_args->GetNumDescriptions(),
                               fmt_args->GetNumAlignments(),
                               *scope,
                               opt.GetMatrixName(),
                               fmt_args->ShowGis(),
                               fmt_args->DisplayHtmlOutput(),
                               opt.GetQueryGeneticCode(),
                               opt.GetDbGeneticCode(),
                               opt.GetSumStatisticsMode(),
                               g_cmdLineArgs->ExecuteRemotely(),
                               db_adapter->GetFilteringAlgorithm(),
                               fmt_args->GetCustomOutputFormatSpec(),
                               g_cmdLineArgs->GetTask() == "megablast",
                               opt.GetMBIndexLoaded());
                
        ///
        /// Use the CLocalBlast class to run a BLAST search 
        ///
        double blast_call_time;
        
        if (g_timingEnabled) {
            blast_call_time = MPI_Wtime();
            gettimeofday(&blastcallStartTime, NULL); /// Wall-clock time
            getrusage(RUSAGE_SELF, &ru_blastcall);   /// Process time
            blastcallStart_u_Time = ru_blastcall.ru_utime;
            blastcallStart_s_Time = ru_blastcall.ru_stime;
            LOG << g_logMsg << "blast call starts,"
                << blast_call_time << ","
                << blastcallStartTime.tv_sec * 1000000
                 + blastcallStartTime.tv_usec << ","
                << blastcallStart_u_Time.tv_sec * 1000000
                 + blastcallStart_u_Time.tv_usec << ","
                << blastcallStart_s_Time.tv_sec * 1000000
                 + blastcallStart_s_Time.tv_usec << ","
                << dbFileName << "," << g_mapCallNo << "," << g_MPI_procName
                << "," << qBlockStart << endl;
            LOG.flush();
        }
        
        /*** Process the input ***/
        for (; !input.End(); formatter.ResetScopeHistory()) {
            CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
            /// GetAllSeqs() instead GetNextSeqBatch()?
            CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
            CLocalBlast lcl_blast(queries, g_optsHndl, db_adapter);
            lcl_blast.SetNumberOfThreads(g_cmdLineArgs->GetNumThreads());
            CRef<CSearchResultSet> results;
            //////////////////////////
            results = lcl_blast.Run();
            //////////////////////////
                        
            /// Blastn original printing
            //ITERATE(CSearchResultSet, result, *results) {
                //formatter.PrintOneResultSet(**result, query_batch);
            //}
            
            for (uint32_t i = 0; i < results->GetNumResults(); ++i) {
                if ((*results)[i].HasAlignments()) {
                    CConstRef<CSeq_align_set> aln_set = (*results)[i].GetSeqAlign();
                    assert(!aln_set->IsEmpty());
                    
                    ITERATE(CSeq_align_set::Tdata, itr_res, aln_set->Get()) {
                        const CSeq_align& s = **itr_res;     
                        
                        ///
                        /// Ref: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/classCSeq__align.html
                        ///
                        /// eScore_BitScore: BLAST-specific bit score
                        /// eScore_EValue: BLAST-specific e-value
                        /// eScore_AlignLength: not a score per se, but a useful metric nonetheless. This is the sum of all aligned segments and all gaps; this excludes introns and discontinuities
                        /// eScore_IdentityCount
                        /// eScore_MismatchCount
                        /// eScore_PercentIdentity_Gapped
                        /// eScore_PercentCoverage
                        /// eScore_PercentIdentity = eScore_PercentIdentity_Gapped 
                        ///                
                        string queryID      = s.GetSeq_id(QUERY).GetSeqIdString();
                        string subID        = s.GetSeq_id(SUBJECT).GetSeqIdString();
                        uint32_t qStart     = s.GetSeqStart(QUERY);
                        uint32_t qEnd       = s.GetSeqStop(QUERY);
                        uint32_t sStart     = s.GetSeqStart(SUBJECT);
                        uint32_t sEnd       = s.GetSeqStop(SUBJECT);
                        qEnd += 1;          /// 0-9 -> 1-10
                        sEnd += 1;          /// 0-9 -> 1-10
                             
                        //uint32_t alignLenGap   = s.GetAlignLength();
                        //uint32_t alignLenUngap = s.GetAlignLength(false);                    
                        
                        ///
                        /// Retrieve def line 
                        ///
                        const CBioseq_Handle& query_bh2 = 
                            scope->GetBioseqHandle(s.GetSeq_id(QUERY));
                        _ASSERT(query_bh2);
                        string defLine = sequence::GetTitle(query_bh2);
                        
                        /// Get GI    
                        vector<string> vecTokens;
                        boost::split(vecTokens, defLine, boost::is_any_of("|"));
                        uint32_t intGi = boost::lexical_cast<uint32_t>(vecTokens[1]);
                        string strGi = vecTokens[1];

                        ///
                        /// Get cutting coords    
                        /// Format: orig_header_chunkID_x_y_v_w
                        ///
                        ///          lower    upper   lower case
                        ///         |xxxxxxx|XXXXXXX|xxxxxxx|
                        ///         x       v       w       y
                        ///
                        ///         x: cutStart,   y: cutEnd
                        ///         v: upperStart, w: upperEnd       
                        /// 
                        ///  chunkID Type-0:     XXX
                        ///          Type-2,4,6: xxxXXX
                        ///          Type-3,5,7: xxxXXXxxx
                        ///          Type-1:     XXXxxx
                        ///
                        string coord = vecTokens[vecTokens.size()-1];
                        assert(coord.length() > 0);
                        vecTokens.clear();
                        boost::split(vecTokens, coord, boost::is_any_of("_"));
                        assert(vecTokens.size() >= 8);
                        int numTokens = vecTokens.size();
                        
                        ///
                        /// read header from mason simulator
                        /// origDefline_readId_pairId(0/1)_strand(f/r)_readLen_
                        /// cutStrat_cutEnd_readStart_readEnd
                        /// ex) >seq1_87_1_r_100_5601_5702_5601_5701
                        ///
                        string cReadId = vecTokens[numTokens-8];
                        string cPairId = vecTokens[numTokens-7];
                        string cStrand = vecTokens[numTokens-6];
                        //uint32_t cLen           
                            //= boost::lexical_cast<uint32_t>(vecTokens[numTokens-5]);
                        uint32_t cutStart        
                            = boost::lexical_cast<uint32_t>(vecTokens[numTokens-4]); /// original cut start (starting from 0)
                        uint32_t cutEnd          
                            = boost::lexical_cast<uint32_t>(vecTokens[numTokens-3]); /// original cut end
                        uint32_t upperStart      
                            = boost::lexical_cast<uint32_t>(vecTokens[numTokens-2]); /// simulated read start (starting from 0)
                        uint32_t upperEnd        
                            = boost::lexical_cast<uint32_t>(vecTokens[numTokens-1]); /// simulated read end (could be longer or shorter than "original cut end"
                         
                        /// Just in case, remove out the hits in the lowercase parts
                        if ((qStart + cutStart < upperStart && qEnd + cutStart < upperStart) ||
                            (qStart + cutStart > upperEnd && qEnd + cutStart > upperEnd)) {
                            if (g_bLogEnabled) {
                                uint32_t cutEnd = boost::lexical_cast<uint32_t>(vecTokens[vecTokens.size()-3]);
                                LOG << g_logMsg << "Warning: A HSP is found in the flank areas!" << endl; 
                                LOG << g_logMsg << "gi, sid, qstart, qend, sstart, send, cutstart, cutend, upperstart, upperend\n";
                                LOG << g_logMsg 
                                    << intGi << ","
                                    << subID << ","
                                    << qStart << ","
                                    << qEnd << ","
                                    << sStart << ","
                                    << sEnd << ","
                                    << cutStart << ","
                                    << cutEnd << ","
                                    << upperStart << ","
                                    << upperEnd << endl;
                                LOG.flush();
                            }
                        }    
                        else {
                            ///
                            /// score, bit_score, evalue...
                            ///
                            int nIdentBases = -1; /// number of identical bases
                            int score = 0, nSum = 0;
                            double bitScore = 0.0, evalue = 0.0;
                            list<int> useThisGi;
                            CAlignFormatUtil::GetAlnScores(s, score, bitScore, evalue, 
                                                           nSum, nIdentBases, useThisGi);
                                                                    
                            ///
                            /// Convert Std-seg and Dense-diag alignments to Dense-seg.
                            /// Ref: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/tabular_8cpp-source.html
                            ///      http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/classCSeq__id.html
                            ///
                            bool isQueryNucl = false, isSubjectNucl = false;
                            const CBioseq_Handle& bioSeqQueryHndl = 
                                scope->GetBioseqHandle(s.GetSeq_id(QUERY));
                            _ASSERT(bioSeqQueryHndl);
                            isQueryNucl = bioSeqQueryHndl.IsNa();
                            
                            const CBioseq_Handle& bioSeqSubjectHndl = 
                                scope->GetBioseqHandle(s.GetSeq_id(SUBJECT));
                            _ASSERT(bioSeqSubjectHndl);
                            isSubjectNucl = bioSeqSubjectHndl.IsNa();
 
                            const bool kTranslated = s.GetSegs().IsStd();
                            CRef<CSeq_align> finalAln(0);                    
                            if (kTranslated) {
                                CRef<CSeq_align> densegAln = s.CreateDensegFromStdseg();
                                if (isQueryNucl && isSubjectNucl)
                                    finalAln = densegAln->CreateTranslatedDensegFromNADenseg();
                                else
                                    finalAln = densegAln;
                            }
                            else if (s.GetSegs().IsDendiag()) {
                                finalAln = CAlignFormatUtil::CreateDensegFromDendiag(s);
                            }   
                            const CDense_seg& ds = (finalAln ? finalAln->GetSegs().GetDenseg() :
                                                    s.GetSegs().GetDenseg());
                                                    
                            CRef<CAlnVec> alnVec;
                            if (!kTranslated && ds.IsSetStrands() && 
                                ds.GetStrands().front() == eNa_strand_minus) {
                                CRef<CDense_seg> reversed_ds(new CDense_seg);
                                reversed_ds->Assign(ds);
                                reversed_ds->Reverse();
                                alnVec.Reset(new CAlnVec(*reversed_ds, *scope));   
                            } else {
                                alnVec.Reset(new CAlnVec(ds, *scope));
                            }    
 
                            ///
                            /// Ref: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/score__builder_8cpp-source.html
                            /// pct_identity = 100.0f * float(*identities) / count_aligned;
                            ///
                            //int alignLen = 0, nGaps = 0, nGapOpens = 0;
                            //CAlignFormatUtil::GetAlignLengths(*alnVec, alignLen, 
                                                            //nGaps, nGapOpens);
                            //float origPercIdent = (alignLen > 0 ? 
                                           //((float)nIdentBases)/alignLen * 100 : 0)
                            //int nMismatches = alignLen - nIdentBases - nGaps;
                            
                            ///
                            /// Doug's filtering
                            ///
                            /// "Percent Identity” looks at every position in the aligned sequences 
                            /// and counts the number that have the same base or amino acid.  
                            /// This count is then divided by the aligned length and then multiplied by 100.
                            ///
                            /// identity = # of identical bases / length of read
                            /// coverage = (read end – read begin)/ length of read 
                            ///
                            /// identity > 50%
                            /// coverage > 90%
                            ///
                            /// Ref: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/score__builder_8cpp-source.html
                            /// pct_coverage = 100.0f * float(covered_bases) / float(seq_len);
                            ///
                            string querySeq = "";
                            string subjectSeq = "";
                            alnVec->SetGapChar('-');
                            alnVec->GetWholeAlnSeqString(QUERY, querySeq);
                            alnVec->GetWholeAlnSeqString(SUBJECT, subjectSeq);
                            assert(querySeq.size() > 0 && subjectSeq.size() > 0);
                            
                            /// Convert to the original coord in the query sequence
                            uint32_t qStartOrig = qStart + cutStart;
                            uint32_t qEndOrig = qEnd + cutStart;               
                            uint32_t scanStart = 0;
                                               
                            ///
                            ///     q.start                        q.end
                            ///        |xxxxxxxxxxxxxxxxxxxxxxxxxxxxx|
                            ///
                            ///    |xxxxxxxxxxx|XXXXXXXXXXXXXXX|xxxxxxxxxxx|
                            /// cutStart   upperStart      upperEnd      cutEnd
                            ///
                            ///         ------>|               |<-----
                            ///             scanStart       scanEnd
                            ///
                            
                            /// Why take min()? Ref: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/tabular_8cpp-source.html
                            uint32_t scanEnd = min(querySeq.size(), subjectSeq.size());
                            if (qStartOrig < upperStart) 
                                scanStart += (upperStart - qStartOrig);
                            if (qEndOrig > upperEnd) 
                                scanEnd -= (qEndOrig - upperEnd);
            
                            ///
                            /// Compare aligned query and subject sequences
                            /// Should consider the gaps in the align when counting 
                            /// num_ident_upperpart and the upperStart and upperEnd 
                            /// which represents the uppercase bases in the query
                            ///
                            uint32_t scanStartInAlign = 0;
                            uint32_t charCount = 0;
                            if (scanStart > 0) {
                                while (charCount < (upperStart - qStartOrig)) {
                                    if (querySeq[scanStartInAlign] != '-') ++charCount;
                                    ++scanStartInAlign;
                                }
                            }
                            charCount = 0;
                            uint32_t scanEndInAlign = min(querySeq.size(), subjectSeq.size());
                            if (scanEnd < min(querySeq.size(), subjectSeq.size())) {
                                while (charCount < (qEndOrig - upperEnd)) {
                                    if (querySeq[scanEndInAlign-1] != '-') ++charCount;
                                    --scanEndInAlign;
                                }
                            }
                            uint32_t identUpperpartLen = 0;                    
                            for (uint32_t i = scanStartInAlign; i < scanEndInAlign; ++i) {
                                if (querySeq[i] == subjectSeq[i]) ++identUpperpartLen;
                            }         
                            
                            float percIdent = float(identUpperpartLen) / 
                                                   (upperEnd - upperStart) * 100;
                                   
                            ///                 
                            /// Adjust qstart and qend for computing coverage
                            ///
                            uint32_t newQStart = qStartOrig, newQEnd = qEndOrig;
                            if (qStartOrig < upperStart) newQStart = upperStart;
                            if (qEndOrig > upperEnd) newQEnd = upperEnd;
                            float percCover = float(newQEnd - newQStart) / 
                                                   (upperEnd - upperStart) * 100;                   
 
                            ///
                            /// Add a HSP to kv and emit
                            ///
                            if (percIdent >= g_identCutoff && percCover >= g_coverCutoff) {
                                /// KV key = intGi
                                ///
                                //structBlRes_t res;
                                //res.subjectId     = boost::lexical_cast<uint32_t>(subID);
                                //res.qStart        = qStart;
                                //res.qEnd          = qEnd;
                                //res.sStart        = sStart;
                                //res.sEnd          = sEnd;
                                //res.eValue        = evalue;
                                //res.bitScore      = boost::lexical_cast<float>(bitScore);
                                //res.upperStart    = upperStart;
                                //res.upperEnd      = upperEnd;
                                //res.identity      = percIdent;
                                //res.coverage      = percCover;
                                
                                //uint32_t newKey   = intGi;
                                //kv->add((char*)&newKey, sizeof(uint32_t), (char*)&res,
                                        //sizeof(structBlRes_t));   
                                    
                                /// KV key = gi_readid_cutX_cutY_F/R_P0/P1
                                ///
                                string strKey = strGi + "_" + 
                                                cReadId + "_" + 
                                                //boost::lexical_cast<string>(upperStart) + "_" +  
                                                //boost::lexical_cast<string>(upperEnd) + "_" + 
                                                boost::lexical_cast<string>(cutStart) + "_" +  /// orginal cutStart from genome seq
                                                boost::lexical_cast<string>(cutEnd) + "_" +    /// orginal cutEnd from genome seq
                                                cStrand + "_" + 
                                                cPairId;
                                
                                structBlResMason_t res;
                                res.subjectId = boost::lexical_cast<uint32_t>(subID);
                                res.qStart    = qStart;
                                res.qEnd      = qEnd;   /// remember this value was increased by one for following Python ( ] range
                                res.sStart    = sStart; 
                                res.sEnd      = sEnd;   /// remember this value was increased by one for following Python ( ] range
                                res.eValue    = evalue;
                                res.bitScore  = boost::lexical_cast<float>(bitScore);
                                res.identity  = percIdent;
                                res.coverage  = percCover;                                

                                const char* newKey = (char*)((strKey).c_str()); 
                                kv->add((char*)newKey, strlen(newKey) + 1, (char*)&res, 
                                        sizeof(structBlResMason_t));
                            }
                        }  
                    }  
                }                       
            }  
        }  
    } CATCH_ALL(status)
    
    assert(status == BLAST_EXIT_SUCCESS);
    
    ///
    /// BLAST call end time
    ///
    if (g_timingEnabled) {
        double blast_call_etime = MPI_Wtime();
        gettimeofday(&blastcallEndTime, NULL); /// Wall-clock time
        getrusage(RUSAGE_SELF, &ru_blastcall); /// process time
        blastcallEnd_u_Time = ru_blastcall.ru_utime;
        blastcallEnd_s_Time = ru_blastcall.ru_stime;
        LOG << g_logMsg << "blast call ends,"
            << blast_call_etime << ","
            << blastcallEndTime.tv_sec * 1000000
             + blastcallEndTime.tv_usec << ","
            << blastcallEnd_u_Time.tv_sec * 1000000
             + blastcallEnd_u_Time.tv_usec << ","
            << blastcallEnd_s_Time.tv_sec * 1000000
             + blastcallEnd_s_Time.tv_usec << ","
            << dbFileName << "," << g_mapCallNo << "," << g_MPI_procName
            << "," << qBlockStart << endl;
        LOG.flush();
    }
}

 

/** Sort function - Passed to MR-MPI sort_values() for sorting blast result
 * string by bit score.
 * @param e1
 * @param e2
 */

inline bool compare_evalue(structEValue_t e1,
                           structEValue_t e2)
{
    return (e1.subjectId != e2.subjectId) ? (e1.subjectId < e2.subjectId) :
                ((e1.eValue != e2.eValue) ? (e1.eValue < e2.eValue) : 
                (e1.bitScore < e2.bitScore));
}

inline void mr_reduce_save_hits( char *key,
                                 int keybytes,
                                 char *multivalue,
                                 int nvalues,
                                 int *valuebytes,
                                 KeyValue *kv,
                                 void *ptr )
{
    ///
    /// Write to file in binary format
    ///
    ofstream outputBinFile((g_hitFileName+".bin").c_str(), ios::binary | ios::app);
    for (size_t n = 0; n < (unsigned)nvalues; ++n) {
        structBlResMason_t* res = (structBlResMason_t*)multivalue;
        ///
        /// KV key = gi_readid_cutX_cutY_F/R_P0/P1
        ///
        vector<string> vecTokens;
        string strKey = key;
        boost::split(vecTokens, strKey, boost::is_any_of("_"));
        assert(vecTokens.size() == 6);
        //cout << key << endl;
        structBlResToSaveHitsMason_t hit;
        hit.gi         = boost::lexical_cast<uint32_t>(vecTokens[0]);
        hit.readId     = boost::lexical_cast<uint32_t>(vecTokens[1]); 
        hit.readStrand = boost::lexical_cast<char>(vecTokens[4]);     /// 'f' or 'r'
        hit.readPairId = boost::lexical_cast<char>(vecTokens[5]);	  /// '0' or '1'
        hit.subjectId  = res->subjectId;
        hit.qStart     = res->qStart;
        hit.qEnd       = res->qEnd;
        hit.sStart     = res->sStart;
        hit.sEnd       = res->sEnd;
        hit.eValue     = res->eValue;
        hit.bitScore   = res->bitScore;
        hit.upperStart = boost::lexical_cast<uint32_t>(vecTokens[2]);
        hit.upperEnd   = boost::lexical_cast<uint32_t>(vecTokens[3]);
        hit.identity   = res->identity;
        hit.coverage   = res->coverage;
        outputBinFile.write((char*)&hit, sizeof(hit));
        multivalue += sizeof(structBlResMason_t);
    }
    outputBinFile.close();    
}

/** Sort by evalue - Passed to MR-MPI reduce() for sorting KMVs by evalue.
 * @param key
 * @param keybytes
 * @param multivalue - collected blast result strings.  
 * @param nvalues
 * @param valuebytes
 * @param kv
 * @param ptr
 */

inline void mr_reduce_sort_multivalues_by_evalue(char *key,
                                   int keybytes,
                                   char *multivalue,
                                   int nvalues,
                                   int *valuebytes,
                                   KeyValue *kv,
                                   void *ptr)
{
    ///
    /// Make structEValue_t = {structBlRes_t* pRec; float evalue;}
    /// and sort by evalue
    ///
    vector<structEValue_t> vecHit;
    for (size_t n = 0; n < nvalues; ++n) {
        structBlRes_t* res = (structBlRes_t*)multivalue;
        structEValue_t strctEvalue;
        strctEvalue.pRec = res;
        strctEvalue.subjectId = res->subjectId;
        strctEvalue.eValue = res->eValue;
        strctEvalue.bitScore = res->bitScore;
        vecHit.push_back(strctEvalue);
        multivalue += sizeof(structBlRes_t);
    }
    
    ///////////////////////////////////////////////////
    sort(vecHit.begin(), vecHit.end(), compare_evalue);
    ///////////////////////////////////////////////////
    
    ///
    /// Write to file in binary format
    ///
    ofstream outputBinFile((g_hitFileName+".bin").c_str(), ios::binary | ios::app);
    for (size_t n = 0; n < (unsigned)nvalues; ++n) {
        structBlRes_t* res = (structBlRes_t*)(vecHit[n].pRec);
        structBlResToSaveHits_t hit;
        hit.gi            = *(uint32_t*)key;
        hit.subjectId      = res->subjectId;
        hit.qStart         = res->qStart;
        hit.qEnd           = res->qEnd;
        hit.sStart         = res->sStart;
        hit.sEnd           = res->sEnd;
        hit.eValue         = res->eValue;
        hit.bitScore       = res->bitScore;
        hit.upperStart     = res->upperStart;
        hit.upperEnd       = res->upperEnd;
        hit.identity       = res->identity;
        hit.coverage       = res->coverage;
        outputBinFile.write((char*)&hit, sizeof(hit));
    }
    outputBinFile.close();
    vecHit.clear();
}


/** Collect MPI node names and ranks
 * @param rank
 * @param nProcs
 * @param mpiComm
 */
 
inline void mpi_collect_node_name(int rank, int nProcs, MPI_Comm mpiComm)
{
    MPI_Status MPI_status;
    char procName[MAXSTR];
    int rankNo;
    int tag = 12345;
    if (rank != MPI_UNDEFINED) { /// if not -32766
        if (rank == 0) {
            for (int p = 1; p < nProcs; ++p) {
                MPI_Recv(&rankNo, 1, MPI_INT, p, tag, mpiComm, &MPI_status);
                MPI_Recv(&procName, MAXSTR, MPI_CHAR, p, tag, mpiComm, &MPI_status);
                g_multimapProcNameRank.insert(pair<string, int>(procName, rankNo));
                g_mapRankProcName.insert(pair<int, string>(rankNo, procName));
            }
        }
        else {
            MPI_Send(&g_MPI_worldRank, 1, MPI_INT, 0, tag, mpiComm);
            MPI_Send(&g_MPI_procName, MAXSTR, MPI_CHAR, 0, tag, mpiComm);
        }
    }
}



/// EOF

