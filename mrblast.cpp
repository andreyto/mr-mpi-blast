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
//  Last updated: 12/19/2011
//
////////////////////////////////////////////////////////////////////////////////

#include "mrblast.hpp"


class CMrMpiBlastApplication : public CNcbiApplication
{
private:
    virtual void Init(void);
    virtual int  Run(void); 
};

void CMrMpiBlastApplication::Init(void)
{    
    if (g_bIsProtein) 
        g_cmdLineArgs.Reset(new CBlastpAppArgs());
    else 
        g_cmdLineArgs.Reset(new CBlastnAppArgs());
        
    HideStdArgs(fHideLogfile | fHideConffile | fHideFullVersion | 
                fHideXmlHelp | fHideDryRun);
                
    SetupArgDescriptions(g_cmdLineArgs->SetCommandLine());  
}

int CMrMpiBlastApplication::Run(void) 
{     
    /// Allow the fasta reader to complain on invalid sequence input
    SetDiagPostLevel(eDiag_Warning);
    
    const CArgs& args = GetArgs();
    RecoverSearchStrategy(args, g_cmdLineArgs);
    CRef<CBlastOptionsHandle> optsHndl(&*g_cmdLineArgs->SetOptions(args));       
    optsHndl->Validate();
    g_optsHndl = optsHndl;
    
    ///
    /// Collect MPI node names and rank number for mapstyle=3 scheduler
    ///
    mpi_collect_node_name(g_MPI_worldRank, g_MPI_nProcs, MPI_COMM_WORLD);
   
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
    ifstream dbListFile(g_dbListFileName.c_str(), ios::in);
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
        boost::split(tok, line, boost::is_any_of("\t"));
        assert(tok.size() == 3);        
        structQueryIndex_t queryIndex;
        queryIndex.qStart      = boost::lexical_cast<uint64_t>(tok[0]);
        queryIndex.qLength     = boost::lexical_cast<uint32_t>(tok[1]);
        queryIndex.uniqQueryId = boost::lexical_cast<uint64_t>(tok[2]);
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
            structBlockBeginLoc_t blockLoc;
            blockLoc.blockBeginLoc = g_vecQueryIndex[q].qStart;
            if (g_bIsQidGi) /// if gi is used for uniq qid
                blockLoc.qIdStart = g_vecQueryIndex[q].uniqQueryId;
            else 
                blockLoc.qIdStart = q+1;
            g_vecBlockBeginLoc.push_back(blockLoc);
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
            aWorkItem.dbNo       = d;
            aWorkItem.blockBegin = g_vecBlockBeginLoc[b].blockBeginLoc;
            aWorkItem.blockEnd   = g_vecBlockBeginLoc[b + 1].blockBeginLoc - 1;
            aWorkItem.qIdStart   = g_vecBlockBeginLoc[b].qIdStart;
            g_vecWorkItem.push_back(aWorkItem);
        }        
        structWorkItem_t aWorkItem;
        aWorkItem.dbNo       = d;
        aWorkItem.blockBegin = g_vecBlockBeginLoc[g_nQueryBlocks - 1].blockBeginLoc;
        aWorkItem.blockEnd   = g_realFileSize;
        aWorkItem.qIdStart   = g_vecBlockBeginLoc[g_nQueryBlocks - 1].qIdStart;
        g_vecWorkItem.push_back(aWorkItem);
    }
    g_nWorkItems = g_vecWorkItem.size();
    
    if (g_MPI_worldRank == 0) {
        cout << "Query file name = " << g_queryFileName << endl;
        cout << "Query index file name = " << g_indexFileName << endl;
        cout << "Database name = " << g_dbName << endl;
        cout << "Database list file name = " << g_dbListFileName << endl;
        cout << "Query block size (bp) = " << g_blockSize << endl;
        cout << "Number of query blocks = " << g_nQueryBlocks << endl;
        cout << "Number of database fragments = " << g_nDbFiles << endl;
        cout << "Number of total work items = " << g_nWorkItems << endl;
        cout << "Number of iterations = " << g_nIter << endl;
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
    if (nRemains) 
        ++nSubWorkItemSets; 

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
    
    ///
    /// Clean up
    ///
    g_vecDbFile.clear();
    g_vecQueryIndex.clear();
    g_vecBlockBeginLoc.clear();
    g_vecWorkItem.clear();
    
    if (!g_searchDatabase.IsNull()) 
        g_searchDatabase.Release();
    if (g_bLogEnabled || g_timingEnabled) 
        g_logFileStream.close();
    if (g_memmapQueryFile.is_open()) 
        g_memmapQueryFile.close(); /// close mmapped file  
    

    return 0;
}
 
 
int main(int argc, char** argv)
{
    ///
    /// Read conf file, mrblast.ini and set parameters
    ///
    ifstream config(CONF_FILE_NAME.c_str(), ios::in);
    if (!config) {
        MPI_Abort(MPI_COMM_WORLD, 1);
        cout << "ERROR: failed to open mr-mpi-blast configuration file, mrblast.ini" << endl;
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
            g_mapStyle       = boost::lexical_cast<int>(parameters["MAPSTYLE"]);
            
            g_bLogEnabled    = boost::lexical_cast<bool>(parameters["LOGENABLED"]);
            g_timingEnabled  = boost::lexical_cast<int>(parameters["TIMING"]);
            g_logFileName    = parameters["LOGFNAME"];
            g_optDumpEnabled = boost::lexical_cast<int>(parameters["OPTDUMP"]);
            
            g_blockSize      = boost::lexical_cast<uint32_t>(parameters["BLOCKSIZE"]);
            g_nIter          = boost::lexical_cast<int>(parameters["NUMITER"]);
            g_bIsProtein     = boost::lexical_cast<bool>(parameters["ISPROTEIN"]);
            g_bIsQidGi       = boost::lexical_cast<bool>(parameters["ISQIDGI"]);
                        
            g_queryFileName  = parameters["QUERYFILENAME"];
            g_indexFileName  = parameters["INDEXFILENAME"];
            g_dbName         = parameters["DATABASENAME"];
            g_dbListFileName = parameters["DBLISTFILENAME"];
            g_outFilePrefix  = parameters["OUTFILEPREFIX"];            
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
    
    ///
    /// MPI setup
    ///
    int MPI_procNameLen;    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_MPI_worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &g_MPI_nProcs);
    MPI_Get_processor_name(g_MPI_procName, &MPI_procNameLen);
    
    ///
    /// Set the "-db -" parameter by default. It's dummy. The actual list of 
    /// databases is read from a database list file, like "dblist.txt".
    /// Check if user sets "-dbsize" if not get the effective db size by using
    /// "blastdbcmd" and set it using "-dbsize".
    ///
    int i;    
    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "-db") == 0) {
            MPI_Abort(MPI_COMM_WORLD, 1);
            cout << "ERROR: '-db' option is not compatible with mr-mpi-blast.\n";
            exit(1);
        }
    }
    char *defArg1 = "-db";
    char *defVal1 = "-";    
    bool bIsDbsizeSet = false;
    /// Check if user sets "-dbsize"
    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "-dbsize") == 0) {
            bIsDbsizeSet = true;
            break;
        }
    }
    
    /// Duplicate the orig argv
    char **newArgv = (char **) malloc((argc+5) * sizeof (char *));
    for (i = 0; i < argc; i++) {
        int len = strlen(argv[i]);
        newArgv[i] = (char *) malloc (len + 1);
        strcpy(newArgv[i], argv[i]); 
    }
    
    /// Add '-db' and '-'
    newArgv[argc] = (char *) malloc (strlen(defArg1) + 1);
    strcpy(newArgv[argc], defArg1);
    newArgv[argc++][strlen(defArg1)+1] = '\0';
    newArgv[argc] = (char *) malloc (strlen(defVal1) + 1);
    strcpy(newArgv[argc], defVal1);
    newArgv[argc++][strlen(defVal1)+1] = '\0';
    newArgv[argc] = NULL;
    
    /// if "dbsize" is not set by user, get the effective
    /// db size of the database and set the number by "-dbsize" BLAST option.
    if (!bIsDbsizeSet) {
        /// Get the effective db size of "g_dbName"
        CRef<CSeqDBExpert> BlastDb;
        CSeqDB::ESeqType seqtype;
        seqtype = g_bIsProtein ? CSeqDB::eProtein : CSeqDB::eNucleotide;;
        BlastDb.Reset(new CSeqDBExpert(g_dbName, seqtype));
        string strDbLen = NStr::UInt8ToString(BlastDb->GetTotalLength());        
        
        /// Bcast the effective dbsize
        char dbLen[50];
        strcpy(dbLen, strDbLen.c_str());
        dbLen[strDbLen.length()+1] = '\0';
        /// Just send dbLen as MPI_CHAR instead MPI_INT_SOMETHING
        MPI_Bcast((char*)dbLen, strDbLen.length()+1, MPI_CHAR, 0, MPI_COMM_WORLD);

        /// Add "-dbsize" option and value
        char *defArg2 = "-dbsize";
        newArgv[argc] = (char *) malloc (strlen(defArg2) + 1);
        strcpy(newArgv[argc], defArg2);
        newArgv[argc++][strlen(defArg2)+1] = '\0';
        newArgv[argc] = (char *) malloc (strlen(dbLen) + 1);
        strcpy(newArgv[argc], dbLen);
        newArgv[argc++][strlen(dbLen)+1] = '\0';
        newArgv[argc] = NULL;
    }
        
    /// Execute main application function
    int ret = CMrMpiBlastApplication().AppMain(argc, newArgv);    
    
    MPI_Finalize();   
    
    return ret;
}


/** mrmpi_blast - MapReduce fuction for BLAST search
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
    pMr->keyalign  = sizeof(uint64_t);  /// The key is query GI
    pMr->mapstyle  = g_mapStyle;        /// master/slave=2, custom scheduler=3
    pMr->outofcore = g_outOfCore;       /// "-1" to disable out-of-core operation
    MPI_Barrier(MPI_COMM_WORLD);

    for (size_t iter = 0; iter < nSubWorkItemSets; ++iter) {
        if (g_MPI_worldRank == 0) 
            cout << "Iteration: " << iter << endl;
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
        pMr->reduce(&mr_reduce_sort_and_save_generic, NULL);
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


/** mr_map_run_blast - MR-MPI Map function - settup NCBI C++ Toolkit env and call blast
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
           
    int dbno;
    uint64_t qBlockStart, qBlockEnd, qIdStart;    
    if (g_nIter > 1) {
        dbno        = g_vecWorkItem[itask + nWorkItemsPerIter * iter].dbNo;
        qBlockStart = g_vecWorkItem[itask + nWorkItemsPerIter * iter].blockBegin;
        qBlockEnd   = g_vecWorkItem[itask + nWorkItemsPerIter * iter].blockEnd;
        qIdStart    = g_vecWorkItem[itask + nWorkItemsPerIter * iter].qIdStart;
    }
    else {
        dbno = g_vecWorkItem[itask].dbNo;
        qBlockStart = g_vecWorkItem[itask].blockBegin;
        qBlockEnd   = g_vecWorkItem[itask].blockEnd;
        qIdStart    = g_vecWorkItem[itask].qIdStart;
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

    ///
    /// Load query from memmapped file
    ///
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
            if (!g_searchDatabase.IsNull()) 
                g_searchDatabase.Release();
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
        if (scope.Empty()) 
            scope.Reset(new CScope(*CObjectManager::GetInstance()));
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
                        string queryID  = s.GetSeq_id(QUERY).GetSeqIdString();
                        string subID    = s.GetSeq_id(SUBJECT).GetSeqIdString();
                        bool bSubIdIsGi = s.GetSeq_id(SUBJECT).IsGi();
                        
                        /// This loc values start from 0~
                        /// NOTE: to make printed hit results those should be 
                        /// started from 1~
                        uint32_t qStart = s.GetSeqStart(QUERY);
                        uint32_t qEnd   = s.GetSeqStop(QUERY);
                        uint32_t sStart = s.GetSeqStart(SUBJECT);
                        uint32_t sEnd   = s.GetSeqStop(SUBJECT);
                       
                        /// if both strands are not matched, sstart and send should change 
                        /// there printing order in the final result 
                        uint8_t qStrand = s.GetSeqStrand(QUERY);
                        uint8_t sStrand = s.GetSeqStrand(SUBJECT);
                        
                        ///
                        /// Get align len with/wo gap(s)
                        ///
                        //uint32_t alignLenGap   = s.GetAlignLength();
                        //uint32_t alignLenUngap = s.GetAlignLength(false);                    
                         
                        ///
                        /// Retrieve query and subject def lines
                        ///
                        //const CBioseq_Handle& query_bh2 = 
                            //scope->GetBioseqHandle(s.GetSeq_id(QUERY));
                        //string qDefLine = sequence::GetTitle(query_bh2);
                        //vector<string> vecTokens;
                        //boost::split(vecTokens, qDefLine, boost::is_any_of("|"));
                        //uint64_t intGi = boost::lexical_cast<uint64_t>(vecTokens[1]);
                        //string strGi = vecTokens[1];
                        
                        string subDefLine;
                        string subjectIdForPrint;
                        if (!bSubIdIsGi) {
                            /// if subject id is not gi, use the defline as subid
                            const CBioseq_Handle& subject_bh = 
                                scope->GetBioseqHandle(s.GetSeq_id(SUBJECT));
                            _ASSERT(subject_bh);
                            subDefLine = sequence::GetTitle(subject_bh);
                            
                            vector<string> vecTokens;
                            boost::split(vecTokens, subDefLine, boost::is_any_of(" "));
                            subjectIdForPrint = vecTokens[0];
                        }
                        else subjectIdForPrint = subID;
 
                        ///
                        /// score, bit_score, evalue...
                        ///
                        int nIdentBases = -1; /// number of identical bases
                        int score = 0, nSum = 0;
                        double bitScore = 0.0;
                        double evalue = 0.0;
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
                        int alignLen = 0, nGaps = 0, nGapOpens = 0;
                        CAlignFormatUtil::GetAlignLengths(*alnVec, alignLen, nGaps, nGapOpens);
                        float origPercIdent = (alignLen > 0 ? 
                            ((float)nIdentBases)/alignLen * 100 : 0);
                        int nMismatches = alignLen - nIdentBases - nGaps;
                                                
                        ///
                        /// Add HSPs to kv and emit
                        ///
 
                        /// Fields to print: 
                        /// query id, subject id, % identity, alignment length, 
                        /// mismatches, gap opens, q. start, q. end, s. start, s. end, 
                        /// evalue, bit score
                        ///
                        uint64_t qIdForPrint;
                        if (!g_bIsQidGi) 
                            qIdForPrint = qIdStart+boost::lexical_cast<uint64_t>(queryID)-1;
                        else
                            qIdForPrint = qIdStart;
                    
                        /// For DEBUG
                        //cout << qIdForPrint << "\t" 
                             //<< subjectIdForPrint << "\t" 
                             //<< origPercIdent << "\t" 
                             //<< alignLen << "\t"
                             //<< nMismatches << "\t"
                             //<< nGaps << "\t"
                             //<< qStart+1 << "\t"
                             //<< qEnd+1 << "\t";                        
                        //if (qStrand != sStrand)
                            //cout << sEnd+1 << "\t" << sStart+1 << "\t";
                        //else 
                            //cout << sStart+1 << "\t" << sEnd+1 << "\t";
                        //cout << evalue << "\t" 
                             //<< bitScore << "\t" 
                             //<< endl;
                        ///
                        
                        structBlResGeneric_t res;
                        strncpy(res.subjectId, subjectIdForPrint.c_str(), MAXSTR-1);
                        res.subjectId[MAXSTR] = '\0';
                        res.identity    = origPercIdent;
                        res.alignLen    = alignLen;
                        res.nMismatches = nMismatches;
                        res.nGaps       = nGaps;
                        res.qStart      = qStart;
                        res.qEnd        = qEnd;                                  
                        if (qStrand != sStrand) {
                            res.sStart  = sEnd+1; 
                            res.sEnd    = sStart+1;   
                        }
                        else {
                            res.sStart  = sStart+1; 
                            res.sEnd    = sEnd+1;   
                        }
                        res.eValue      = evalue;
                        res.bitScore    = bitScore;
                        
                        const char* newKey = (char*)(&qIdForPrint); 
                        kv->add((char*)newKey, sizeof(uint64_t), (char*)&res, 
                                sizeof(structBlResGeneric_t)); 
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

 
/** mr_reduce_sort_and_save_generic - Sort by evalue and save hits in binary.
 * Passed to MR-MPI reduce() for sorting KMVs by evalue.
 * @param key
 * @param keybytes
 * @param multivalue - collected blast result strings.  
 * @param nvalues
 * @param valuebytes
 * @param kv
 * @param ptr
 */

inline void mr_reduce_sort_and_save_generic(char *key,
                                            int keybytes,
                                            char *multivalue,
                                            int nvalues,
                                            int *valuebytes,
                                            KeyValue *kv,
                                            void *ptr)
{
    ///
    /// Make structEValue_t = {structBlResGeneric_t* pRec; float evalue; ...}
    /// and sort by evalue
    ///
    vector<structEValue_t> vecHit;
    for (size_t n = 0; n < nvalues; ++n) {
        structBlResGeneric_t* res = (structBlResGeneric_t*)multivalue;
        structEValue_t strctEvalue;
        strctEvalue.pRec = res;
        strncpy(strctEvalue.subjectId, res->subjectId, MAXSTR);
        strctEvalue.eValue = res->eValue;
        strctEvalue.bitScore = res->bitScore;
        vecHit.push_back(strctEvalue);
        multivalue += sizeof(structBlResGeneric_t);
    }
    
    ///////////////////////////////////////////////////
    sort(vecHit.begin(), vecHit.end(), compare_evalue_generic);
    ///////////////////////////////////////////////////
    
    ///
    /// Write to file in binary format
    ///
    ofstream outputBinFile((g_hitFileName+".bin").c_str(), ios::binary | ios::app);
    for (size_t n = 0; n < (unsigned)nvalues; ++n) {
        structBlResGeneric_t* res = (structBlResGeneric_t*)(vecHit[n].pRec);
        structBlResToSaveHits_t hit;
        hit.queryId     = *(uint64_t*)key;
        strncpy(hit.subjectId, res->subjectId, MAXSTR);
        hit.identity    = res->identity;
        hit.alignLen    = res->alignLen;
        hit.nMismatches = res->nMismatches;
        hit.nGaps       = res->nGaps;
        hit.qStart      = res->qStart;
        hit.qEnd        = res->qEnd;
        hit.sStart      = res->sStart;
        hit.sEnd        = res->sEnd;
        hit.eValue      = res->eValue;
        hit.bitScore    = res->bitScore;
        outputBinFile.write((char*)&hit, sizeof(hit));
        
        /// DEBUG
        //cout << *(uint64_t*)key << "\t" 
             //<< res->subjectId << "\t" 
             //<< res->identity << "\t"              
             //<< res->alignLen << "\t"
             //<< res->nMismatches << "\t"
             //<< res->nGaps << "\t"
             //<< res->qStart << "\t"
             //<< res->qEnd << "\t";             
             //<< res->sStart << "\t" 
             //<< res->sEnd << "\t";
             //<< res->eValue << "\t" 
             //<< res->bitScore << "\t" 
             //<< endl;
        ///
    }
    outputBinFile.close();
    vecHit.clear();
}

/** compare_evalue_generic - Sort function - Passed to MR-MPI sort_values() for sorting blast result
 * string by bit score.
 * @param e1
 * @param e2
 */

inline bool compare_evalue_generic(structEValue_t e1,
                                   structEValue_t e2)
{
    ///
    /// Possible sort criteria
    /// Expect Value: default
    /// Max Score: By the bit score of HSPs, similar to Expect Value
    /// Query Coverage: By the percent of length coverge for the query
    /// Max Identity: By the maximal percent ID of the HSPs
    ///
    int ret = strcmp(e1.subjectId, e2.subjectId);
    return (ret != 0) ? ((e1.eValue != e2.eValue) ? (e1.eValue < e2.eValue) : 
                (e1.bitScore < e2.bitScore)) : (ret);
}


/** mpi_collect_node_name - Collect MPI node names and ranks
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

