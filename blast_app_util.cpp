/*  $Id: blast_app_util.cpp 188930 2010-04-15 20:23:30Z maning $
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
 * Author: Christiam Camacho
 *
 */

/** @file blast_app_util.cpp
 *  Utility functions for BLAST command line applications
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: blast_app_util.cpp 188930 2010-04-15 20:23:30Z maning $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <ncbi_pch.hpp>
#include "blast_app_util.hpp"

#include <serial/serial.hpp>
#include <serial/objostr.hpp>

#include <objtools/data_loaders/blastdb/bdbloader.hpp>
#include <algo/blast/api/remote_blast.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>     // for CObjMgr_QueryFactory
#include <algo/blast/api/blast_options_builder.hpp>
#include <algo/blast/api/search_strategy.hpp>
#include <algo/blast/blastinput/blast_input.hpp>    // for CInputException
#include <algo/blast/blastinput/psiblast_args.hpp>
#include <algo/blast/blastinput/tblastn_args.hpp>
#include <algo/blast/blastinput/blast_scope_src.hpp>
#include <objmgr/util/sequence.hpp>

BEGIN_NCBI_SCOPE
USING_SCOPE(objects);
USING_SCOPE(blast);

CRef<blast::CRemoteBlast> 
InitializeRemoteBlast(CRef<blast::IQueryFactory> queries,
                      CRef<blast::CBlastDatabaseArgs> db_args,
                      CRef<blast::CBlastOptionsHandle> opts_hndl,
                      bool verbose_output,
                      const string& client_id /* = kEmptyStr */,
                      CRef<objects::CPssmWithParameters> pssm 
                        /* = CRef<objects::CPssmWithParameters>() */)
{
    _ASSERT(queries || pssm);
    _ASSERT(db_args);
    _ASSERT(opts_hndl);

    CRef<CRemoteBlast> retval;

    CRef<CSearchDatabase> search_db = db_args->GetSearchDatabase();
    if (search_db.NotEmpty()) {
        if (pssm.NotEmpty()) {
            _ASSERT(queries.Empty());
            retval.Reset(new CRemoteBlast(pssm, opts_hndl, *search_db));
        } else {
            retval.Reset(new CRemoteBlast(queries, opts_hndl, *search_db));
        }
    } else {
        if (pssm.NotEmpty()) {
            NCBI_THROW(CInputException, eInvalidInput,
                       "Remote PSI-BL2SEQ is not supported");
        } else {
            // N.B.: there is NO scope needed in the GetSubjects call because
            // the subjects (if any) should have already been added in 
            // InitializeSubject 
            retval.Reset(new CRemoteBlast(queries, opts_hndl,
                                         db_args->GetSubjects()));
        }
    }
    if (verbose_output) {
        retval->SetVerbose();
    }
    if (client_id != kEmptyStr) {
        retval->SetClientId(client_id);
    }
    return retval;
}

void
InitializeSubject(CRef<blast::CBlastDatabaseArgs> db_args, 
                  CRef<blast::CBlastOptionsHandle> opts_hndl,
                  bool is_remote_search,
                  CRef<blast::CLocalDbAdapter>& db_adapter, 
                  CRef<objects::CScope>& scope)
{
    db_adapter.Reset();

    _ASSERT(db_args.NotEmpty());
    CRef<CSearchDatabase> search_db = db_args->GetSearchDatabase();

    // Initialize the scope... 
    if (is_remote_search) {
        const bool is_protein = 
            Blast_SubjectIsProtein(opts_hndl->GetOptions().GetProgramType())
            ? true : false;
        SDataLoaderConfig config(is_protein);
        if (search_db.NotEmpty()) {
            config.m_BlastDbName = search_db->GetDatabaseName();
        }
        CBlastScopeSource scope_src(config);
        // configure scope to fetch sequences remotely for formatting
        if (scope.NotEmpty()) {
            scope_src.AddDataLoaders(scope);
        } else {
            scope = scope_src.NewScope();
        }
    } else {
        if (scope.Empty()) {
            scope.Reset(new CScope(*CObjectManager::GetInstance()));
        }
    }
    _ASSERT(scope.NotEmpty());

    // ... and then the subjects
    CRef<IQueryFactory> subjects;
    if ( (subjects = db_args->GetSubjects(scope)) ) {
        _ASSERT(search_db.Empty());
        db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl));
    } else {
        _ASSERT(search_db.NotEmpty());
        try { 
            // Try to open the BLAST database even for remote searches, as if
            // it is available locally, it will be better to fetch the
            // sequence data for formatting from this (local) source
            CRef<CSeqDB> seqdb = search_db->GetSeqDb();
            db_adapter.Reset(new CLocalDbAdapter(*search_db));
            scope->AddDataLoader(RegisterOMDataLoader(seqdb));
        } catch (const CSeqDBException&) {
            // The BLAST database couldn't be found, report this for local
            // searches, but for remote searches go on.
            if (is_remote_search ) {
                db_adapter.Reset(new CLocalDbAdapter(*search_db));
            } else {
                throw;
            }
        }
    }
}

string RegisterOMDataLoader(CRef<CSeqDB> db_handle)
{
    // the blast formatter requires that the database coexist in
    // the same scope with the query sequences
    CRef<CObjectManager> om = CObjectManager::GetInstance();
    CBlastDbDataLoader::RegisterInObjectManager(*om, db_handle, true,
                        CObjectManager::eDefault,
                        CBlastDatabaseArgs::kSubjectsDataLoaderPriority);
    CBlastDbDataLoader::SBlastDbParam param(db_handle);
    return CBlastDbDataLoader::GetLoaderNameFromArgs(param);
}

/// Real implementation of search strategy extraction
/// @todo refactor this code so that it can be reused in other contexts
static void
s_ExportSearchStrategy(CNcbiOstream* out,
                     CRef<blast::IQueryFactory> queries,
                     CRef<blast::CBlastOptionsHandle> options_handle,
                     CRef<blast::CBlastDatabaseArgs> db_args,
                     CRef<objects::CPssmWithParameters> pssm 
                       /* = CRef<objects::CPssmWithParameters>() */)
{
    if ( !out ) {
        return;
    }
    _ASSERT(db_args);
    _ASSERT(options_handle);

    try { 
        CRef<CRemoteBlast> rmt_blast =
            InitializeRemoteBlast(queries, db_args, options_handle, false,
                                  kEmptyStr, pssm);
        CRef<CBlast4_request> req = rmt_blast->GetSearchStrategy(); 
        // N.B.: If writing XML, be sure to call SetEnforcedStdXml on the
        // stream!
        *out << MSerial_AsnText << *req;
    } catch (const CBlastException& e) {
        if (e.GetErrCode() == CBlastException::eNotSupported) {
            NCBI_THROW(CInputException, eInvalidInput, 
                       "Saving search strategies with gi lists is currently "
                       "not supported");
        }
        throw;
    }
}

/// Converts a list of Bioseqs into a TSeqLocVector. All Bioseqs are added to
/// the same CScope object
/// @param subjects Bioseqs to convert
static TSeqLocVector
s_ConvertBioseqs2TSeqLocVector(const CBlast4_subject::TSequences& subjects)
{
    TSeqLocVector retval;
    CRef<CScope> subj_scope(new CScope(*CObjectManager::GetInstance()));
    ITERATE(CBlast4_subject::TSequences, bioseq, subjects) {
        subj_scope->AddBioseq(**bioseq);
        CRef<CSeq_id> seqid = FindBestChoice((*bioseq)->GetId(),
                                             CSeq_id::BestRank);
        const TSeqPos length = (*bioseq)->GetInst().GetLength();
        CRef<CSeq_loc> sl(new CSeq_loc(*seqid, 0, length-1));
        retval.push_back(SSeqLoc(sl, subj_scope));
    }
    return retval;
}

/// Import PSSM into the command line arguments object
static void 
s_ImportPssm(const CBlast4_queries& queries,
             CRef<blast::CBlastOptionsHandle> opts_hndl,
             blast::CBlastAppArgs* cmdline_args)
{
    CRef<CPssmWithParameters> pssm
        (const_cast<CPssmWithParameters*>(&queries.GetPssm()));
    CPsiBlastAppArgs* psi_args = NULL;
    CTblastnAppArgs* tbn_args = NULL;

    if ( (psi_args = dynamic_cast<CPsiBlastAppArgs*>(cmdline_args)) ) {
        psi_args->SetInputPssm(pssm);
    } else if ( (tbn_args = 
                 dynamic_cast<CTblastnAppArgs*>(cmdline_args))) {
        tbn_args->SetInputPssm(pssm);
    } else {
        EBlastProgramType p = opts_hndl->GetOptions().GetProgramType();
        string msg("PSSM found in saved strategy, but not supported ");
        msg += "for " + Blast_ProgramNameFromType(p);
        NCBI_THROW(CBlastException, eNotSupported, msg);
    }
}

/// Import queries into command line arguments object
static void 
s_ImportQueries(const CBlast4_queries& queries,
                CRef<blast::CBlastOptionsHandle> opts_hndl,
                blast::CBlastAppArgs* cmdline_args)
{
    CRef<CTmpFile> tmpfile(new CTmpFile(CTmpFile::eNoRemove));

    // Stuff the query bioseq or seqloc list in the input stream of the
    // cmdline_args
    if (queries.IsSeq_loc_list()) {
        const CBlast4_queries::TSeq_loc_list& seqlocs =
            queries.GetSeq_loc_list();
        CFastaOstream out(tmpfile->AsOutputFile(CTmpFile::eIfExists_Throw));
        out.SetFlag(CFastaOstream::eAssembleParts);
        
        EBlastProgramType prog = opts_hndl->GetOptions().GetProgramType();
        SDataLoaderConfig dlconfig(!!Blast_QueryIsProtein(prog));
        dlconfig.OptimizeForWholeLargeSequenceRetrieval();
        CBlastScopeSource scope_src(dlconfig);
        CRef<CScope> scope(scope_src.NewScope());

        ITERATE(CBlast4_queries::TSeq_loc_list, itr, seqlocs) {
            CBioseq_Handle bh = scope->GetBioseqHandle(**itr);
            CConstRef<CBioseq> bioseq = bh.GetCompleteBioseq();
            out.Write(*bioseq);
        }
        scope.Reset();
        scope_src.RevokeBlastDbDataLoader();

    } else {
        _ASSERT(queries.IsBioseq_set());
        const CBlast4_queries::TBioseq_set& bioseqs =
            queries.GetBioseq_set();
        CFastaOstream out(tmpfile->AsOutputFile(CTmpFile::eIfExists_Throw));
        out.SetFlag(CFastaOstream::eAssembleParts);

        ITERATE(CBioseq_set::TSeq_set, seq_entry, bioseqs.GetSeq_set()){
            out.Write(**seq_entry);
        }
    }

    const string& fname = tmpfile->GetFileName();
    tmpfile.Reset(new CTmpFile(fname));
    cmdline_args->SetInputStream(tmpfile);
}

/// Import the database and return it in a CBlastDatabaseArgs object
static CRef<blast::CBlastDatabaseArgs>
s_ImportDatabase(const CBlast4_subject& subj, 
                 CBlastOptionsBuilder& opts_builder,
                 bool subject_is_protein,
                 bool is_remote_search)
{
    _ASSERT(subj.IsDatabase());
    CRef<CBlastDatabaseArgs> db_args(new CBlastDatabaseArgs);
    const CSearchDatabase::EMoleculeType mol = subject_is_protein
        ? CSearchDatabase::eBlastDbIsProtein
        : CSearchDatabase::eBlastDbIsNucleotide;
    const string dbname(subj.GetDatabase());
    CRef<CSearchDatabase> search_db(new CSearchDatabase(dbname, mol));

    if (opts_builder.HaveEntrezQuery()) {
        string limit(opts_builder.GetEntrezQuery());
        search_db->SetEntrezQueryLimitation(limit);
        if ( !is_remote_search ) {
            string msg("Entrez query '");
            msg += limit + string("' will not be processed locally.\n");
            msg += string("Please use the -remote option.");
            throw runtime_error(msg);
        }
    }

    if (opts_builder.HaveGiList()) {
        CSeqDBGiList *gilist = new CSeqDBGiList();
        ITERATE(list<int>, gi, opts_builder.GetGiList()) {
            gilist->AddGi(*gi);
        }
        search_db->SetGiList(gilist);
    }

    if (opts_builder.HasDbFilteringAlgorithmId()) {
        int algo_id = opts_builder.GetDbFilteringAlgorithmId();
        search_db->SetFilteringAlgorithm(algo_id);
    }

    db_args->SetSearchDatabase(search_db);
    return db_args;
}

/// Import the subject sequences into a CBlastDatabaseArgs object
static CRef<blast::CBlastDatabaseArgs>
s_ImportSubjects(const CBlast4_subject& subj, bool subject_is_protein)
{
    _ASSERT(subj.IsSequences());
    CRef<CBlastDatabaseArgs> db_args(new CBlastDatabaseArgs);
    TSeqLocVector subjects = 
        s_ConvertBioseqs2TSeqLocVector(subj.GetSequences());
    CRef<CScope> subj_scope = subjects.front().scope;
    CRef<IQueryFactory> subject_factory(new CObjMgr_QueryFactory(subjects));
    db_args->SetSubjects(subject_factory, subj_scope, subject_is_protein);
    return db_args;
}

/// Imports search strategy, using CImportStrategy.
static void
s_ImportSearchStrategy(CNcbiIstream* in, 
                       blast::CBlastAppArgs* cmdline_args,
                       bool is_remote_search, 
                       bool override_query, 
                       bool override_subject)
{
    if ( !in ) {
        return;
    }

    CRef<CBlast4_request> b4req;
    try { 
        b4req = ExtractBlast4Request(*in);
    } catch (const CSerialException&) {
        NCBI_THROW(CInputException, eInvalidInput, 
                   "Failed to read search strategy");
    }

    CImportStrategy strategy(b4req);

    CRef<blast::CBlastOptionsHandle> opts_hndl = strategy.GetOptionsHandle();
    cmdline_args->SetOptionsHandle(opts_hndl);
    const EBlastProgramType prog = opts_hndl->GetOptions().GetProgramType();
    cmdline_args->SetTask(strategy.GetTask());

    // Get the subject
    if (override_subject) {
        ERR_POST(Warning << "Overriding database/subject in saved strategy");
    } else {
        CRef<blast::CBlastDatabaseArgs> db_args;
        CRef<CBlast4_subject> subj = strategy.GetSubject();
    const bool subject_is_protein = Blast_SubjectIsProtein(prog) ? true : false;

        if (subj->IsDatabase()) {
            CBlastOptionsBuilder bob(strategy.GetProgram(), strategy.GetService(), CBlastOptions::eBoth);
            bob.GetSearchOptions(&strategy.GetAlgoOptions(), &strategy.GetProgramOptions());
            db_args = s_ImportDatabase(*subj, bob, subject_is_protein,
                                       is_remote_search);
        } else {
            db_args = s_ImportSubjects(*subj, subject_is_protein);
        }
        _ASSERT(db_args.NotEmpty());
        cmdline_args->SetBlastDatabaseArgs(db_args);
    }

    // Get the query, queries, or pssm
    if (override_query) {
        ERR_POST(Warning << "Overriding query in saved strategy");
    } else {
        CRef<CBlast4_queries> queries = strategy.GetQueries();
        if (queries->IsPssm()) {
            s_ImportPssm(*queries, opts_hndl, cmdline_args);
        } else {
            s_ImportQueries(*queries, opts_hndl, cmdline_args);
        }
        // Set the range restriction for the query, if applicable
        const TSeqRange query_range = strategy.GetQueryRange();
        if (query_range != TSeqRange::GetEmpty()) {
            cmdline_args->GetQueryOptionsArgs()->SetRange(query_range);
        }
    }

}

void
RecoverSearchStrategy(const CArgs& args, blast::CBlastAppArgs* cmdline_args)
{
    CNcbiIstream* in = cmdline_args->GetImportSearchStrategyStream(args);
    if ( !in ) {
        return;
    }
    const bool is_remote_search = 
        (args[kArgRemote].HasValue() && args[kArgRemote].AsBoolean());
    const bool override_query = (args[kArgQuery].HasValue() && 
                                 args[kArgQuery].AsString() != kDfltArgQuery);
    const bool override_subject = CBlastDatabaseArgs::HasBeenSet(args);
    s_ImportSearchStrategy(in, cmdline_args, is_remote_search, override_query,
                           override_subject);
    if (CMbIndexArgs::HasBeenSet(args)) {
        ERR_POST(Warning << "Overriding megablast BLAST DB indexed options in saved strategy");
    }
}

// Process search strategies
// FIXME: save program options,
// Save task if provided, no other options (only those in the cmd line) should
// be saved
void
SaveSearchStrategy(const CArgs& args,
                   blast::CBlastAppArgs* cmdline_args,
                   CRef<blast::IQueryFactory> queries,
                   CRef<blast::CBlastOptionsHandle> opts_hndl,
                   CRef<objects::CPssmWithParameters> pssm 
                     /* = CRef<objects::CPssmWithParameters>() */)
{
    CNcbiOstream* out = cmdline_args->GetExportSearchStrategyStream(args);
    if ( !out ) {
        return;
    }

    s_ExportSearchStrategy(out, queries, opts_hndl, 
                           cmdline_args->GetBlastDatabaseArgs(), 
                           pssm);
}

END_NCBI_SCOPE
