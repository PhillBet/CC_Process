
from Bio.Seq import Seq
from libs import utils
import libs.columns as CL


def csv_format(data):

    reg_entity = {
        '_id': data.get('regulated_data', None).get('regulated_id', None),
        'type': data.get('regulated_data', None).get('regulated_type', None)
    }
    # print(reg_entity)
    first_gene = utils.get_distance_to_first_gene(
        site_id=data.get('site_data', None).get('site_id', None),
        reg_entity=reg_entity,
        regulated_genes=utils.regulated_genes(reg_entity)
    )
    # print(first_gene)

    ri_seq = data.get('site_data', None).get('site_sequence', None)
    if ri_seq and data.get('regulated_data', None).get('regulated_strand', None):
        if data.get('regulated_data', None).get('regulated_strand', None) == 'reverse':
            ri_seq = Seq(ri_seq)
            ri_seq = ri_seq.reverse_complement()

    row_dict = {
        CL.RI_CYC_ID: data.get('ri_cyc_id', None),
        CL.RI_ID: data.get('_id', None),
        CL.TF_CYC_ID: data.get('tf_data', None).get('tf_cyc_id', None),
        CL.TF_ID: data.get('tf_data', None).get('tf_id', None),
        CL.TF_NAME: data.get('tf_data', None).get('tf_name', None),
        CL.FINAL_STATE: data.get('regulator_data', None).get('regulator_name', None),
        CL.REG_ENTITY_CYC_ID: data.get('regulated_data', None).get('regulated_cyc_id', None),
        CL.REG_ENTITY_ID: data.get('regulated_data', None).get('regulated_id', None),
        CL.REG_ENTITY_NAME: data.get('regulated_data', None).get('regulated_name', None),
        CL.REG_ENTITY_STRAND: data.get('regulated_data', None).get('regulated_strand', None),
        CL.REG_ENTITY_TYPE: data.get('regulated_data', None).get('regulated_type', None),
        CL.SITE_CYC_ID: data.get('site_cyc_id', None),
        CL.SITE_ID: data.get('site_data', None).get('site_id', None),
        CL.SITE_LEFT: data.get('site_data', None).get('site_left', None),
        CL.SITE_RIGHT: data.get('site_data', None).get('site_right', None),
        CL.SITE_LENGTH: data.get('site_data', None).get('site_length', None),
        CL.RI_FUNCTION: data.get('function', None),
        CL.RI_DIST_TO_TSS: data.get('relativeDistSitePromoter', None),
        CL.RI_DIST_FIRST_GENE: first_gene.get('distance'),
        CL.RI_FIRST_GENE_ID: first_gene.get('_id'),
        CL.RI_ORIENTATION: data.get('orientation', None),
        CL.RI_SEQUENCE: ri_seq,
    }
    return row_dict
