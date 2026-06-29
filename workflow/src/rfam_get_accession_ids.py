import argparse
import json
import csv
import pandas
import Bio.Entrez as Entrez

RFAM_FULL_REGION_COLS = ["family", "seq_accession", "start", "end", "bitscore", "evalue", "cm_start", "cm_end", "truncated", "type"]

import time
from Bio import Entrez

def map_nuc_to_assembly(accession_list):
    """
    Maps a list of nucleotide accessions to their parent genome assembly accessions.
    """
    mapping_results = {}
    
    # Convert list of accessions into a comma-separated string
    id_string = ",".join(accession_list)
    
    try:
        # Link from the nucleotide database ('nuccore') to the assembly database ('assembly')
        with Entrez.elink(dbfrom="nuccore", db="assembly", id=id_string, idtype="acc") as handle:
            record = Entrez.read(handle)
            
        # Parse the nested XML response
        for entry in record:
            # Extract the original query accession
            query_acc = entry["IdList"][0] 
            mapping_results[query_acc] = []
            
            # Check if any linked assembly IDs were found
            if entry["LinkSetDb"]:
                for link in entry["LinkSetDb"][0]["Link"]:
                    assembly_id = link["Id"]
                    
                    # Convert the internal NCBI UID to a standard GCF/GCA accession number
                    with Entrez.esummary(db="assembly", id=assembly_id) as sum_handle:
                        sum_record = Entrez.read(sum_handle)
                        assembly_acc = sum_record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
                        mapping_results[query_acc].append(assembly_acc)
            else:
                mapping_results[query_acc].append("No linked assembly found")
                
    except Exception as e:
        print(f"An error occurred: {e}")
        
    return mapping_results

def seq_accession_to_gid(seq_accessions):
    search_term = " OR ".join(seq_accessions)
    search_handle = Entrez.esearch(db="nucleotide", term=search_term)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    gi_ids = search_results.get("IdList", [])
    print(f"Found {len(gi_ids)}/{len(seq_accessions)} UIDs")

    return gi_ids

def gi_id_to_assembly_uid(gi_ids):
    link_handle = Entrez.elink(
        dbfrom="nucleotide",
        db="assembly",
        id=",".join(gi_ids),
        linkname="nuccore_assembly"
    )
    link_results = Entrez.read(link_handle)
    link_handle.close()
    
    assembly_uids = []
    
    for link_set in link_results:
        assembly_uids.append([link['Id'] for link in link_set['LinkSetDb'][0]['Link']])

    return assembly_uids

def seq_accession_to_assembly_id(seq_accessions):
    id_string = ",".join(seq_accessions)
    with Entrez.elink(dbfrom="nuccore", db="assembly", id=id_string, idtype="acc") as handle:
        record = Entrez.read(handle)

    seqacc_to_uid = {}

    # from what I've seen omly one record is returned but maybe if you increase batch size??

    for entry in record:
        query_seq_acc = entry["IdList"]
        uids = [link['Id'] for link in entry['LinkSetDb'][0]['Link']]

        if len(query_seq_acc) != len(uids):
            print(f"WARNING: sequence accessions do not map one to one to uids. {len(query_seq_acc)} query accessions and {len(uids)} genome accession ids found")

        seqacc_to_uid.update(zip(query_seq_acc, uids))

    return seqacc_to_uid


def uid_to_assembly(assembly_uids):
    """
    Args: assembly_uids list of genome assembly uids
    Returns: dict mapping asembly uids to accession name (GCA)"""
    summary_handle = Entrez.esummary(db="assembly", id=",".join(assembly_uids))
    summary_results = Entrez.read(summary_handle)
    summary_handle.close()
    
    uid_to_gca = {}
    doc_summaries = summary_results.get("DocumentSummarySet", {}).get("DocumentSummary", [])
    
    if not isinstance(doc_summaries, list):
        doc_summaries = [doc_summaries]
    
    uid_to_gca = {}

    for doc in doc_summaries:
        uid = doc.attributes.get("uid")
        
        gca_acc = doc.get("AssemblyAccession")
        if uid and gca_acc:
            uid_to_gca[uid] = gca_acc

    return uid_to_gca

def get_genome_assembly_acc(sequence_table, chunk_size=200):
    sequence_table = sequence_table.copy()
    assembly_acc = []

    for chunk_idx in range(0, len(sequence_table), chunk_size):
        batch = sequence_table["seq_accession"].iloc[chunk_idx : chunk_idx + chunk_size].tolist()

        seq_to_uid = seq_accession_to_assembly_id(batch)

        uid_to_gca = uid_to_assembly(seq_to_uid.values())

        # for seq in seq_to_uid.keys():
        #     uid = seq_to_uid[seq]
        #     gca = uid_to_gca[uid]

            #sequence_table.loc[sequence_table["seq_accession"] == seq, "genome_accession"] = gca
        #UNFORTUNATELY MAPPING DOESN'T WORK
        assembly_acc.extend(uid_to_gca.values())

    return assembly_acc


def batch_map_sequences_to_assemblies(sequence_table, chunk_size=200):
    """
    Map nucleotide sequence accessions to genome assembly accessions (GCA/GCF) using NCBI Entrez.
    
    Approach:
    1. esearch: text accessions -> internal UIDs
    2. elink: nucleotide UIDs -> assembly UIDs
    3. esummary (assembly): assembly UIDs -> GCA/GCF accessions
    4. esummary (nucleotide): nucleotide UIDs -> original sequence IDs (validation)
    
    Arguments:
        sequence_table: Pandas dataframe with seq_accession column
        chunk_size: Number of sequences per batch (default 200)
    
    Returns:
        Pandas dataframe with added genome_accession column
    """
    
    sequence_table = sequence_table.copy()
    sequence_table["genome_accession"] = None
    
    total_sequences = len(sequence_table)
    
    print(f"\n{'='*80}")
    print(f"Starting batch mapping: {total_sequences} sequences in chunks of {chunk_size}")
    print(f"{'='*80}\n")
    
    for chunk_idx in range(0, total_sequences, chunk_size):
        chunk_start = chunk_idx
        chunk_end = min(chunk_idx + chunk_size, total_sequences)
        chunk_num = (chunk_idx // chunk_size) + 1
        
        chunk_df = sequence_table.iloc[chunk_start:chunk_end]
        chunk_accessions = chunk_df["seq_accession"].tolist()
        
        print(f"\n[Chunk {chunk_num}] Processing indices {chunk_start}-{chunk_end-1} ({len(chunk_accessions)} sequences)")
        print(f"  Accessions: {chunk_accessions[:5]}{'...' if len(chunk_accessions) > 5 else ''}")
        
        #try:
            # ===== STEP A: Convert text accessions to internal UIDs =====
        print(f"  Step A: esearch (text → UIDs)...", end=" ", flush=True)

        search_term = " OR ".join(chunk_accessions)
        search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=chunk_size)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        gi_ids = search_results.get("IdList", [])
        print(f"Found {len(gi_ids)}/{len(chunk_accessions)} UIDs")

        seq_to_gi_ids = dict(zip(chunk_accessions, gi_ids))
        
        if not gi_ids:
            print(f"    ⚠ WARNING: No UIDs found for this chunk")
            time.sleep(0.5 if not Entrez.api_key else 0.2)
            continue
        
        if len(gi_ids) < len(chunk_accessions):
            missing_count = len(chunk_accessions) - len(gi_ids)
            print(f"    ⚠ WARNING: {missing_count} sequences not found in esearch")
        
        # ===== STEP B: Link nucleotide UIDs to assembly UIDs =====
        # STILL CAN'T FIGURE OUT HOW TO GET RELIABLE MAPPING FROM INPUT IDS TO OUTPUT ASSAMBLY IDS
        # I HAVE TRIED DOING THESE ELINK REQUESTS INDIVIDUALLY BUT THE CONNECTION TIMES OUT
        print(f"  Step B: elink (nucleotide → assembly)...", end=" ", flush=True)
        
        link_handle = Entrez.elink(
            dbfrom="nucleotide",
            db="assembly",
            id=",".join(gi_ids),
            linkname="nuccore_assembly"
        )
        link_results = Entrez.read(link_handle)
        link_handle.close()
        
        gi_to_assembly = {}
        assembly_uids = []
        
        for link_set in link_results:
            # Extract the GI IDs and assembly IDs
            gi_ids = link_set['IdList']
            assembly_uids = [link['Id'] for link in link_set['LinkSetDb'][0]['Link']]

            # Create the mapping
            gi_to_assembly.update(dict(zip(gi_ids, assembly_uids)))
        
        print(f"Linked {len(gi_to_assembly.values())/len(gi_ids)} to assemblies")
        
        if not assembly_uids:
            print(f"    ⚠ WARNING: No assembly UIDs found for this chunk")
            time.sleep(0.5 if not Entrez.api_key else 0.2)
            continue
        
        # ===== STEP C: Fetch assembly summaries (UIDs → GCA/GCF) =====
        print(f"  Step C: esummary assembly (UIDs → GCA/GCF)...", end=" ", flush=True)
        summary_handle = Entrez.esummary(db="assembly", id=",".join(assembly_uids))
        summary_results = Entrez.read(summary_handle)
        summary_handle.close()
        
        uid_to_gca = {}
        doc_summaries = summary_results.get("DocumentSummarySet", {}).get("DocumentSummary", [])
        
        if not isinstance(doc_summaries, list):
            doc_summaries = [doc_summaries]
        
        gca_found = 0
        for doc in doc_summaries:
            try:
                # Handle both attribute styles depending on Entrez.read() version
                uid = doc.attributes.get("uid")
                
                gca_acc = doc.get("AssemblyAccession")
                if uid and gca_acc:
                    uid_to_gca[uid] = gca_acc
                    gca_found += 1
            except (KeyError, AttributeError) as e:
                print(f"    ⚠ Error parsing assembly summary: {e}")
                continue
        
        print(f"Resolved {gca_found}/{len(assembly_uids)} UIDs to GCA/GCF")

        print(uid_to_gca)

        for seq_acc, gi_id in seq_to_gi_ids.items():
            if gi_id not in gi_to_assembly:
                print(f"    ⚠ WARNING: GI ID {gi_id} has no assembly link")
                continue
            
            uid = gi_to_assembly[gi_id]
            if uid not in uid_to_gca:
                print(f"    ⚠ WARNING: UID {uid} has no GCA mapping")
                continue
            
            gca = uid_to_gca[uid]
            sequence_table.loc[sequence_table["seq_accession"] == seq_acc, "genome_accession"] = gca
        
        # ===== STEP D: Fetch nucleotide summaries (validation + original accessions) =====
        # print(f"  Step D: esummary nucleotide (validation)...", end=" ", flush=True)
        # nucl_summary_handle = Entrez.esummary(db="nucleotide", id=",".join(gi_ids))
        # nucl_summaries = Entrez.read(nucl_summary_handle)
        # nucl_summary_handle.close()

        # print(nucl_summaries)
        
        # if not isinstance(nucl_summaries, list):
        #     nucl_summaries = [nucl_summaries]
        
        # # ===== STEP E: Build final mapping =====
        # mappings_complete = 0
        # for doc in nucl_summaries:
        #     try:
        #         orig_seq_id = doc.get("Caption")
        #         internal_gi = doc.get("attributes", {}).get("uid") or doc.get("uid")
                
        #         if not orig_seq_id or not internal_gi:
        #             continue
                
        #         linked_assembly_uid = gi_to_assembly.get(internal_gi)
        #         final_gca = uid_to_gca.get(linked_assembly_uid) if linked_assembly_uid else None
                
        #         if orig_seq_id and final_gca:
        #             mapping_results[orig_seq_id] = final_gca
        #             mappings_complete += 1
        #         elif orig_seq_id and not final_gca:
        #             mapping_results[orig_seq_id] = None
            
        #     except (KeyError, AttributeError, TypeError) as e:
        #         print(f"    ⚠ Error processing nucleotide summary: {e}")
        #         continue
        
        # print(f"Complete mappings: {mappings_complete}/{len(gi_ids)}")
            
        # except Exception as e:
        #     print(f"\n  ERROR in chunk {chunk_num}: {e}")
        #     sequence_table.loc[chunk_df.index, "map_status"] = "failed"
        #     raise e
        
        # finally:
        #     # Polite rate limiting
        #     time.sleep(0.5 if not Entrez.api_key else 0.2)
    
    # # ===== Apply results to dataframe =====
    # print(f"\n{'='*80}")
    # print("Applying results to dataframe...")
    
    # for seq_id, gca_accession in mapping_results.items():
    #     # Find rows matching this sequence ID
    #     matching_rows = sequence_table[sequence_table["seq_accession"] == seq_id].index
        
    #     if len(matching_rows) > 0:
    #         for idx in matching_rows:
    #             if gca_accession:
    #                 sequence_table.loc[idx, "genome_accession"] = gca_accession
    #                 sequence_table.loc[idx, "map_status"] = "found"
    #             else:
    #                 sequence_table.loc[idx, "map_status"] = "no_link"
    
    # # Summary statistics
    # found_count = (sequence_table["map_status"] == "found").sum()
    # no_link_count = (sequence_table["map_status"] == "no_link").sum()
    # failed_count = (sequence_table["map_status"] == "failed").sum()
    # pending_count = (sequence_table["map_status"] == "pending").sum()
    
    # print(f"\nResults:")
    # print(f"  ✓ Found mappings:     {found_count}")
    # print(f"  ⚠ No assembly links:  {no_link_count}")
    # print(f"  ✗ Failed:             {failed_count}")
    # print(f"  ○ Pending:            {pending_count}")
    # print(f"  {'─'*40}")
    # print(f"  Total:                {len(sequence_table)}")
    # print(f"{'='*80}\n")
    
    return sequence_table

def map_sequences_to_assemblies(sequence_table, chunk_size=200):
    """
    Map nucleotide sequence accessions to genome assembly accessions (GCA/GCF) using NCBI Entrez.
    
    Approach:
    1. esearch: text accessions → gi_ids (list)
    2. elink: gi_ids → assembly_uids (list)
    3. esummary (assembly): assembly_uids → GCA/GCF accessions
    4. esummary (nucleotide): gi_ids → original sequence IDs + trace chain to GCA
    
    Arguments:
        sequence_table: Pandas dataframe with seq_accession column
        chunk_size: Number of sequences per batch (default 200)
    
    Returns:
        Pandas dataframe with added genome_accession column
    """
    
    sequence_table = sequence_table.copy()
    sequence_table["genome_accession"] = None
    sequence_table["map_status"] = "pending"
    
    total_sequences = len(sequence_table)
    
    print(f"\n{'='*80}")
    print(f"Starting batch mapping: {total_sequences} sequences in chunks of {chunk_size}")
    print(f"{'='*80}\n")
    
    for chunk_idx in range(0, total_sequences, chunk_size):
        chunk_start = chunk_idx
        chunk_end = min(chunk_idx + chunk_size, total_sequences)
        chunk_num = (chunk_idx // chunk_size) + 1
        
        chunk_df = sequence_table.iloc[chunk_start:chunk_end]
        chunk_accessions = chunk_df["seq_accession"].tolist()
        
        print(f"\n[Chunk {chunk_num}] Processing indices {chunk_start}-{chunk_end-1} ({len(chunk_accessions)} sequences)")
        print(f"  Accessions: {chunk_accessions[:5]}{'...' if len(chunk_accessions) > 5 else ''}")
        
        try:
            # ===== STEP A: Convert text accessions to internal UIDs =====
            print(f"  Step A: esearch (text → gi_ids)...", end=" ", flush=True)
            search_term = " OR ".join(chunk_accessions)
            search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=chunk_size)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            gi_ids = search_results.get("IdList", [])
            print(f"Found {len(gi_ids)}/{len(chunk_accessions)} gi_ids")
            
            if not gi_ids:
                print(f"    ⚠ WARNING: No gi_ids found for this chunk")
                sequence_table.loc[chunk_df.index, "map_status"] = "failed"
                time.sleep(0.5 if not Entrez.api_key else 0.2)
                continue
            
            if len(gi_ids) < len(chunk_accessions):
                missing_count = len(chunk_accessions) - len(gi_ids)
                print(f"    ⚠ WARNING: {missing_count} sequences not found in esearch")
            
            # ===== STEP B: Link nucleotide gi_ids to assembly_uids =====
            print(f"  Step B: elink (gi_ids → assembly_uids)...", end=" ", flush=True)
            
            link_handle = Entrez.elink(
                dbfrom="nucleotide",
                db="assembly",
                id=",".join(gi_ids),
                linkname="nuccore_assembly"
            )
            link_results = Entrez.read(link_handle)
            link_handle.close()
            
            assembly_uids = []
            links_found = 0
            
            for link_set in link_results:
                link_set_dbs = link_set.get("LinkSetDb", [])
                
                if link_set_dbs and "Link" in link_set_dbs[0]:
                    try:
                        # Take first (most recent) assembly link
                        assembly_uid = link_set_dbs[0]["Link"][0]["Id"]
                        assembly_uids.append(assembly_uid)
                        links_found += 1
                    except (IndexError, KeyError, TypeError):
                        pass
            
            print(f"Linked {links_found}/{len(gi_ids)} to assembly_uids")
            
            if not assembly_uids:
                print(f"    ⚠ WARNING: No assembly_uids found for this chunk")
                sequence_table.loc[chunk_df.index, "map_status"] = "no_link"
                time.sleep(0.5 if not Entrez.api_key else 0.2)
                continue
            
            # ===== STEP C: Fetch assembly summaries (assembly_uids → GCA/GCF) =====
            print(f"  Step C: esummary assembly (assembly_uids → GCA/GCF)...", end=" ", flush=True)
            summary_handle = Entrez.esummary(db="assembly", id=",".join(assembly_uids))
            summary_results = Entrez.read(summary_handle)
            summary_handle.close()
            
            doc_summaries = summary_results.get("DocumentSummarySet", {}).get("DocumentSummary", [])
            
            if not isinstance(doc_summaries, list):
                doc_summaries = [doc_summaries]
            
            # Build assembly_uid → GCA mapping
            assembly_uid_to_gca = {}
            gca_found = 0
            
            for doc in doc_summaries:
                try:
                    uid = doc.attributes.get("uid")
                    gca_acc = doc.get("AssemblyAccession")
                    
                    if uid and gca_acc:
                        assembly_uid_to_gca[uid] = gca_acc
                        gca_found += 1
                except (KeyError, AttributeError, TypeError) as e:
                    print(f"    ⚠ Error parsing assembly summary: {e}")
                    continue
            
            print(f"Resolved {gca_found}/{len(assembly_uids)} assembly_uids to GCA/GCF")
            
            # ===== STEP D: Fetch nucleotide summaries (gi_ids → original accessions) =====
            print(f"  Step D: esummary nucleotide (validation & tracing)...", end=" ", flush=True)
            nucl_summary_handle = Entrez.esummary(db="nucleotide", id=",".join(gi_ids))
            nucl_summaries = Entrez.read(nucl_summary_handle)
            nucl_summary_handle.close()
            
            if not isinstance(nucl_summaries, list):
                nucl_summaries = [nucl_summaries]
            
            # Build gi_id → original_accession mapping from Caption
            gi_id_to_accession = {}
            for doc in nucl_summaries:
                try:
                    gi_id = doc.attributes.get("uid")
                    orig_accession = doc.get("Caption")
                    
                    if gi_id and orig_accession:
                        gi_id_to_accession[gi_id] = orig_accession
                except (KeyError, AttributeError, TypeError):
                    pass
            
            print(f"Retrieved {len(gi_id_to_accession)}/{len(gi_ids)} original accessions")
            
            # ===== STEP E: Trace chain and apply to dataframe =====
            mappings_complete = 0
            
            for idx, gi_id in enumerate(gi_ids):
                orig_accession = gi_id_to_accession.get(gi_id)
                
                if not orig_accession:
                    print(f"    ⚠ Could not retrieve original accession for gi_id {gi_id}")
                    continue
                
                # Get the assembly UID that was linked to this gi_id
                # Note: elink returns results in same order as input gi_ids
                if idx < len(assembly_uids):
                    assembly_uid = assembly_uids[idx]
                    gca_accession = assembly_uid_to_gca.get(assembly_uid)
                    
                    if gca_accession:
                        sequence_table.loc[
                            sequence_table["seq_accession"] == orig_accession,
                            "genome_accession"
                        ] = gca_accession
                        
                        sequence_table.loc[
                            sequence_table["seq_accession"] == orig_accession,
                            "map_status"
                        ] = "found"
                        
                        mappings_complete += 1
                    else:
                        sequence_table.loc[
                            sequence_table["seq_accession"] == orig_accession,
                            "map_status"
                        ] = "no_gca"
                else:
                    sequence_table.loc[
                        sequence_table["seq_accession"] == orig_accession,
                        "map_status"
                    ] = "no_assembly_link"
            
            print(f"Complete mappings: {mappings_complete}/{len(gi_ids)}")
        
        except Exception as e:
            print(f"\n  ✗ ERROR in chunk {chunk_num}: {e}")
            sequence_table.loc[chunk_df.index, "map_status"] = "failed"
            import traceback
            traceback.print_exc()
        
        finally:
            # Polite rate limiting
            time.sleep(0.5 if not Entrez.api_key else 0.2)
    
    # ===== Summary statistics =====
    print(f"\n{'='*80}")
    print("Results:")
    
    found_count = (sequence_table["map_status"] == "found").sum()
    no_gca_count = (sequence_table["map_status"] == "no_gca").sum()
    no_link_count = (sequence_table["map_status"] == "no_link").sum()
    no_assembly_count = (sequence_table["map_status"] == "no_assembly_link").sum()
    failed_count = (sequence_table["map_status"] == "failed").sum()
    pending_count = (sequence_table["map_status"] == "pending").sum()
    
    print(f"  ✓ Found mappings:        {found_count}")
    print(f"  ⚠ Assembly UID found, no GCA:   {no_gca_count}")
    print(f"  ⚠ No assembly links:     {no_link_count}")
    print(f"  ⚠ No assembly UID:       {no_assembly_count}")
    print(f"  ✗ Failed:                {failed_count}")
    print(f"  ○ Pending:               {pending_count}")
    print(f"  {'─'*40}")
    print(f"  Total:                   {len(sequence_table)}")
    print(f"{'='*80}\n")
    
    return sequence_table

# def get_assembly_from_sequence(seq_id):
#     # Step 1: Find the internal NCBI ID for the nucleotide sequence
#     search_handle = Entrez.esearch(db="nucleotide", term=seq_id)
#     search_results = Entrez.read(search_handle)
#     search_handle.close()
    
#     if not search_results["IdList"]:
#         print(f"No nucleotide record found for: {seq_id}")
#         return None
    
#     gi_id = search_results["IdList"][0]
    
#     # Step 2: Link the nucleotide ID to the assembly database
#     link_handle = Entrez.elink(dbfrom="nucleotide", db="assembly", id=gi_id)
#     link_results = Entrez.read(link_handle)
#     link_handle.close()
    
#     # Extract linked assembly ID
#     link_set_dbs = link_results[0].get("LinkSetDb", [])
#     if not link_set_dbs:
#         print(f"No linked assembly found for sequence ID: {seq_id}")
#         return None
        
#     assembly_id = link_set_dbs[0]["Link"][0]["Id"]
    
#     # Step 3: Fetch the summary of the assembly to get the GCA/GCF accession string
#     summary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
#     summary_results = Entrez.read(summary_handle)
#     summary_handle.close()
    
#     # Extract the textual Accession
#     assembly_acc = summary_results["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
#     return assembly_acc

def get_genome_accession_mapping(sequence_table, db="nucleotide", batch_size=100):
    """Given a list of sequence accession IDs, return a mapping of sequence ID to genome assembly accession ID (GCA/GCF) using Entrez links.
    Arguments:
        sequence_table: Pandas dataframe with seq_accession column
        db: Database to search in (default "nucleotide")
        batch_size: Number of IDs to query in each batch (default 100)
    Returns: Pandas dataframe with additional genome_accession column mapping sequence accession ID to genome assembly accession ID (GCA/GCF)"""
    
    sequence_table = sequence_table.copy()
    sequence_table["genome_accession"] = None
    
    for i in range(0, len(sequence_table), batch_size):
        batch = sequence_table["seq_accession"].iloc[i : i + batch_size].tolist()
        id_string = ",".join(batch)
        
        link_handle = Entrez.elink(
            dbfrom=db, 
            db="assembly", 
            id=id_string, 
            linkname=f"{db}_assembly"
        )
        
        try:
            link_results = Entrez.read(link_handle)
        except Exception as e:
            print(f"error querying entrez: {e}")
            continue
        finally:
            link_handle.close()
        
        # Parse linkset results: each linkset contains one sequence and its linked assemblies
        for linkset in link_results:
            seq_id = linkset.get("IdList", [None])[0]
            assembly_uid = None

            print(seq_id)
            
            # Extract assembly UID from nested LinkSetDb structure
            if seq_id and "LinkSetDb" in linkset:
                for linksetdb in linkset["LinkSetDb"]:
                    if "Link" in linksetdb and len(linksetdb["Link"]) > 0:
                        assembly_uid = linksetdb["Link"][0].get("Id")
                        break
            
            if assembly_uid:
                # Fetch the assembly record to get the actual GCA/GCF accession
                try:
                    summary_handle = Entrez.esummary(db="assembly", id=assembly_uid)
                    summary = Entrez.read(summary_handle)
                    
                    if summary.get("DocumentSummarySet"):
                        assembly_accession = summary["DocumentSummarySet"]["DocumentSummary"][0].get("AssemblyAccession")
                        sequence_table.loc[sequence_table["seq_accession"] == seq_id, "genome_accession"] = assembly_accession
                    else:
                        sequence_table.loc[sequence_table["seq_accession"] == seq_id, "genome_accession"] = None
                except Exception as e:
                    print(f"Error fetching assembly accession for {seq_id}: {e}")
                    sequence_table.loc[sequence_table["seq_accession"] == seq_id, "genome_accession"] = None
                finally:
                    summary_handle.close()
            else:
                # No linked assembly found
                if seq_id:
                    sequence_table.loc[sequence_table["seq_accession"] == seq_id, "genome_accession"] = None
    
    return sequence_table


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input file, rfam full region table")
    parser.add_argument("--family", help="Rfam family")
    parser.add_argument("--output-list", help="Output file, list of accession IDs in JSON format")
    #parser.add_argument("--output-table", help="Output file, metadata table in tsv format")
    parser.add_argument("--entrez-email", help="Email associated with Entrez requests, to comply with policy but seems to be optional")
    parser.add_argument("--entrez-api-key", help="Entrez API key, if provided increases limit from 3 to 10 req/s", required=False)
    args = parser.parse_args()

    Entrez.email = args.entrez_email
    print(f"Setting email: {args.entrez_email}")
    if args.entrez_api_key:
        Entrez.api_key = args.entrez_api_key

    full_regions_tbl = pandas.read_table(args.input, header=None, names=RFAM_FULL_REGION_COLS)

    # filter for the specified family
    family_regions_tbl = full_regions_tbl[full_regions_tbl["family"] == args.family]

    #family_regions_tbl = get_genome_accession_mapping(family_regions_tbl)
    #family_regions_tbl = map_sequences_to_assemblies(family_regions_tbl)
    genome_list = get_genome_assembly_acc(family_regions_tbl)

    # save genome accessions as json list 
    with open(args.output_list, 'w') as f:
        json.dump(genome_list, f, indent=4)

    # save rest of the table
    #family_regions_tbl.to_csv(args.output_table, sep='\t')

if __name__ == "__main__":
    main()


# Example usage:
# python rfam_get_accession_ids.py --input resources/rfam/Rfam.full_region --family RF00001 --output results/00_controls/rfam/RF00001/accession_ids.json