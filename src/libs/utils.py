import csv

import multigenomic_api as mg_api

import libs.columns as CL
import libs.constants as EC
import libs.collections.ri as ri_coll

regulator_types = {
    'product': EC.PD_COLLECTION,
    'regulatoryComplex': EC.RCPLX_COLLECTION,
    'regulatoryContinuant': EC.RCONT_COLLECTION,
}

regulated_types = {
    'promoter': EC.PM_COLLECTION,
    'gene': EC.GN_COLLECTION,
    'transcriptionUnit': EC.TU_COLLECTION,
}


def set_genome_intervals():
    '''
    Set the genes ranges to calculate the closest_genes.

    Param

    Returns
        genes_ranges, List, Array of coordinate pairs of the calculated ranges.
    '''
    genome_length = EC.GENOME_LENGTH
    intervals = EC.INTERVALS
    intervals_length = int(genome_length / intervals)
    genes_ranges = []
    for interval in range(intervals):
        genes_ranges.append(
            [intervals_length + ((interval - 1) * intervals_length), intervals_length + (interval * intervals_length)])
    return genes_ranges


def get_center_pos(left_pos, right_pos):
    '''
    Calculates the center center position of the chromosome.

    Param
        left_pos, String, Start position in the sequence (it's converted to Integer).
        right_pos, String, End position in the sequence (it's converted to Integer).

    Returns
        center_pos, Float, Center position in the sequence.
    '''
    center_pos = int(right_pos) - int(left_pos)
    center_pos = (center_pos / 2) + int(left_pos)
    return center_pos


def regulated_genes(reg_entity):
    regulated_genes = []
    transcription_units = []
    if reg_entity.get('type') == "gene":
        gene = mg_api.genes.find_by_id(reg_entity.get('_id'))
        gene_object = {
            "_id": gene.id,
            "name": gene.name,
        }
        if gene_object not in regulated_genes:
            regulated_genes.append(gene_object)
    if reg_entity.get('type') == "promoter":
        transcription_units = mg_api.transcription_units.find_by_promoter_id(
            reg_entity.get('_id'))
    elif reg_entity.get('type') == "transcriptionUnit":
        trans_unit = mg_api.transcription_units.find_by_id(
            reg_entity.get('_id'))
        transcription_units.append(trans_unit)
    if transcription_units:
        for tu in transcription_units:
            for gene_id in tu.genes_ids:
                gene = mg_api.genes.find_by_id(gene_id)
                gene_object = {
                    "_id": gene.id,
                    "name": gene.name,
                }
                if gene_object not in regulated_genes:
                    regulated_genes.append(gene_object)
    return regulated_genes


def get_first_gene_of_tu(genes, promoter):
    first_gene = None
    if len(genes) > 0:
        gene = genes[0]
        first_gene = mg_api.genes.find_by_id(gene.get("_id"))
        if promoter:
            if promoter.strand == "reverse":
                first_gene.left_end_position = first_gene.right_end_position
            first_gene.left_end_position = first_gene.left_end_position or first_gene.fragments[
                0].left_end_position

            for gene in genes:
                current_gene = mg_api.genes.find_by_id(gene.get("_id"))
                if promoter.strand == "forward":
                    if current_gene.left_end_position:
                        if current_gene.left_end_position < first_gene.left_end_position:
                            first_gene = current_gene
                    elif current_gene.fragments:
                        for fragment in current_gene.fragments:
                            if fragment.left_end_position < first_gene.left_end_position:
                                first_gene = current_gene
                                first_gene.left_end_position = fragment.left_end_position

                elif promoter.strand == "reverse":
                    if current_gene.left_end_position:
                        if current_gene.right_end_position > first_gene.left_end_position:
                            first_gene = current_gene
                            first_gene.left_end_position = first_gene.right_end_position
                    elif current_gene.fragments:
                        for fragment in current_gene.fragments:
                            if fragment.right_end_position > first_gene.left_end_position:
                                first_gene = current_gene
                                first_gene.left_end_position = fragment.right_end_position

                first_gene.left_end_position = first_gene.left_end_position or first_gene.fragments[
                    0].left_end_position
        first_gene = {
            "_id": first_gene.id,
            "leftEndPosition": first_gene.left_end_position
        }
    return first_gene


def get_distance_to_first_gene(site_id, reg_entity, regulated_genes):
    distance_to_first_gene = None
    promoter = None
    first_gene = None
    first_gene_dict = {}
    if site_id and regulated_genes != []:
        reg_sites = mg_api.regulatory_sites.find_by_id(
            site_id)
        # print(reg_entity.get('type'))
        if reg_entity.get('type') == "gene":
            #print(site_id, reg_entity, reg_entity.get('type'), regulated_genes)
            first_gene = mg_api.genes.find_by_id(reg_entity.get('_id'))
            first_gene_id = first_gene.id
            if first_gene.strand:
                if reg_sites.absolute_position:
                    if first_gene.strand == "forward":
                        distance_to_first_gene = reg_sites.absolute_position - first_gene.left_end_position
                    else:
                        distance_to_first_gene = first_gene.right_end_position - reg_sites.absolute_position
                elif reg_sites.left_end_position and reg_sites.right_end_position:
                    absolute_position = (
                        reg_sites.left_end_position + reg_sites.right_end_position) / 2
                    if first_gene.strand == "forward":
                        distance_to_first_gene = absolute_position - \
                            first_gene.left_end_position
                    else:
                        distance_to_first_gene = first_gene.left_end_position - \
                            absolute_position
        else:
            if reg_entity.get('type') == "promoter":
                promoter = mg_api.promoters.find_by_id(
                    reg_entity.get('_id'))
                first_gene = get_first_gene_of_tu(regulated_genes, promoter)
            elif reg_entity.get('type') == "transcriptionUnit":
                trans_unit = mg_api.transcription_units.find_by_id(
                    reg_entity.get('_id'))
                if trans_unit.promoters_id:
                    promoter = mg_api.promoters.find_by_id(
                        trans_unit.promoters_id)
                first_gene = get_first_gene_of_tu(regulated_genes, promoter)
            #print(site_id, reg_entity.get('type'), reg_entity.get('_id'))
            if promoter:
                if reg_sites.absolute_position:
                    if promoter.strand == "forward":
                        distance_to_first_gene = reg_sites.absolute_position - \
                            first_gene["leftEndPosition"]
                    else:
                        distance_to_first_gene = first_gene["leftEndPosition"] - \
                            reg_sites.absolute_position
                elif reg_sites.left_end_position and reg_sites.right_end_position:
                    absolute_position = (
                        reg_sites.left_end_position + reg_sites.right_end_position) / 2
                    if promoter.strand == "forward":
                        distance_to_first_gene = absolute_position - \
                            first_gene["leftEndPosition"]
                    else:
                        distance_to_first_gene = first_gene["leftEndPosition"] - \
                            absolute_position
            first_gene_id = first_gene.get("_id")

        if first_gene:
            first_gene_dict = {
                '_id': first_gene_id,
                'distance': distance_to_first_gene,
            }

    return first_gene_dict


def gen_csv_file(table_name, rows, fields):
    with open(f'{table_name}.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def get_regulated_data(regulated, db, collection_ids):
    regulated_id = regulated.get('_id', None)
    if not regulated_id:
        return {}
    regulated_type = regulated.get('type', None)

    collection_name = regulated_types.get(regulated_type, None)
    # print(regulated_type, collection_name)
    collection = db[collection_name]
    regulated_name = collection.find_one(
        {'_id': regulated_id}).get('abbreviatedName', None)
    if regulated_name is None:
        regulated_name = collection.find_one(
            {'_id': regulated_id}).get('name', None)
    if regulated_name is None:
        regulated_name = regulated.get('name', None)
    if regulated_type == 'transcriptionUnit':
        operon_collection = db['operons']
        operon_id = collection.find_one(
            {'_id': regulated_id}).get('operons_id', None)
        regulated_strand = operon_collection.find_one(
            {'_id': operon_id}).get('strand', None)
    else:
        regulated_strand = collection.find_one(
            {'_id': regulated_id}).get('strand', None)
    regulated_cyc_id = collection_ids.find_one({'_id': regulated_id}).get(
        'objectOriginalSourceId', None)
    regulated_data = {
        'regulated_cyc_id': regulated_cyc_id,
        'regulated_id': regulated_id,
        'regulated_name': regulated_name,
        'regulated_type': regulated_type,
        'regulated_strand': regulated_strand,

    }
    return regulated_data


def get_tf_data(regulator_id, regulator_name, db, collection_ids):
    mg_tf_obj = mg_api.transcription_factors.find_tf_id_by_conformation_id(
        regulator_id)

    if mg_tf_obj:
        # print(mg_tf_obj[0].id)
        tf_cyc_id = collection_ids.find_one({'_id': mg_tf_obj[0].id}).get(
            'objectOriginalSourceId', None)
        tf_data = {
            'tf_cyc_id': tf_cyc_id,
            'tf_id': mg_tf_obj[0].id,
            'tf_name': mg_tf_obj[0].abbreviated_name,
        }
        return tf_data
    else:
        collection_name = 'transcriptionFactors'
        collection = db[collection_name]
        mg_tf_obj = collection.find_one(
            {'abbreviatedName': regulator_name})
        if mg_tf_obj:
            tf_cyc_id = collection_ids.find_one({'_id': mg_tf_obj.get('_id')}).get(
                'objectOriginalSourceId', None)
            tf_data = {
                'tf_cyc_id': tf_cyc_id,
                'tf_id': mg_tf_obj.get('_id'),
                'tf_name': regulator_name,
            }
            return tf_data
        else:
            # print(regulator_name)
            return {}


def get_regulator_data(regulator, db):
    regulator_id = regulator.get('_id', None)
    if not regulator_id:
        return {}
    regulator_type = regulator.get('type', None)
    collection_name = regulator_types.get(regulator_type, None)
    # print(regulator_type, collection_name)
    collection = db[collection_name]
    regulator_name = collection.find_one(
        {'_id': regulator_id}).get('abbreviatedName', None)
    if regulator_name is None:
        regulator_name = collection.find_one(
            {'_id': regulator_id}).get('name', None)
    if regulator_name is None:
        regulator_name = regulator.get('name', None)
    regulator_data = {
        'regulator_id': regulator_id,
        'regulator_type': regulator_type,
        'regulator_name': regulator_name,
    }
    return regulator_data


def get_site_data(site_id, db):
    collection_name = 'regulatorySites'
    collection = db[collection_name]
    if site_id is None:
        return {}
    site_left = collection.find_one(
        {'_id': site_id}).get('leftEndPosition', None)
    site_right = collection.find_one(
        {'_id': site_id}).get('rightEndPosition', None)
    site_length = collection.find_one(
        {'_id': site_id}).get('length', None)
    site_sequence = collection.find_one(
        {'_id': site_id}).get('sequence', None)
    site_data = {
        'site_id': site_id,
        'site_left': site_left,
        'site_right': site_right,
        'site_length': site_length,
        'site_sequence': site_sequence,

    }
    return site_data


def set_csv_format(collection_name, data):
    row_dict = None
    if collection_name == EC.RI_COLLECTION:
        row_dict = ri_coll.csv_format(data)

    return row_dict
