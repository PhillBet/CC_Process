

import pymongo
import multigenomic_api

from libs import utils
from libs import arguments
import libs.columns as CL
import libs.constants as EC


def run(**kwargs):

    multigenomic_api.connect(kwargs.get('mg_database'), kwargs.get('url'))
    mongo_client = pymongo.MongoClient(kwargs.get('url'))
    ids_db = mongo_client[kwargs.get('id_database')]
    db = mongo_client[kwargs.get('mg_database')]

    table_name = kwargs.get('table_name', None)

    if table_name == EC.ALL_TABLES or table_name == EC.RI_TABLE:
        collection_name = EC.RI_COLLECTION
        print(f'Generating {collection_name} collection')
        collection_ids = ids_db['identifiers']
        collection = db[collection_name]
        query = {}
        mg_objects = collection.find(query)
        rows = []
        genes_ranges = utils.set_genome_intervals()
        for obj in mg_objects:  # [0:1000]:
            ri_cyc_id = collection_ids.find_one({'_id': obj.get('_id', None)})
            obj.update({'ri_cyc_id': ri_cyc_id.get(
                'objectOriginalSourceId', None)})
            #print(obj.get('_id', None))

            site_cyc_id = collection_ids.find_one(
                {'_id': obj.get('regulatorySites_id', None)})
            if site_cyc_id is not None:
                obj.update({'site_cyc_id': site_cyc_id.get(
                    'objectOriginalSourceId', None)})

            regulated = obj.get('regulatedEntity', None)
            regulated_data = utils.get_regulated_data(
                regulated, db, collection_ids)
            if regulated_data is not None:
                obj.update({'regulated_data': regulated_data})

            regulator = obj.get('regulator', None)
            regulator_data = utils.get_regulator_data(regulator, db)
            if regulator_data is not None:
                obj.update({'regulator_data': regulator_data})

            tf_data = utils.get_tf_data(
                regulator_data.get('regulator_id', None), regulator_data.get('regulator_name', None), db, collection_ids)
            if tf_data is not None:
                obj.update({'tf_data': tf_data})

            site_data = utils.get_site_data(
                obj.get('regulatorySites_id', None), db)
            if site_data is not None:
                obj.update({'site_data': site_data})
            # print(regulator_data)

            obj.update({'genes_ranges': genes_ranges})

            row = utils.set_csv_format(collection_name, obj)
            if row:
                rows.append(row)
        utils.gen_csv_file(
            table_name=collection_name,
            rows=rows,
            fields=CL.RI_FIELDS
        )
        sites_without_length = []
        for row in rows:  # Sin length o sin positions L/R
            if not row.get(CL.SITE_LENGTH, None) and row.get(CL.SITE_RIGHT) and row.get(CL.SITE_LEFT):
                sites_without_length.append(row)
            if row.get(CL.SITE_LENGTH, None) and not row.get(CL.SITE_RIGHT) and not row.get(CL.SITE_LEFT):
                sites_without_length.append(row)
            if not row.get(CL.SITE_LENGTH, None) and not row.get(CL.SITE_RIGHT) and not row.get(CL.SITE_LEFT):
                sites_without_length.append(row)
        utils.gen_csv_file(
            table_name=f'{collection_name}_no_length',
            rows=sites_without_length,
            fields=CL.RI_FIELDS
        )
        sites_without_tf = []
        for row in rows:  # Sin TF
            if not row.get(CL.TF_CYC_ID, None):
                sites_without_tf.append(row)
        utils.gen_csv_file(
            table_name=f'{collection_name}_no_tf',
            rows=sites_without_tf,
            fields=CL.RI_FIELDS
        )

        sites_without_tf = []
        for row in rows:  # Sin TF
            if not row.get(CL.RI_FIRST_GENE_ID, None):
                sites_without_tf.append(row)
        utils.gen_csv_file(
            table_name=f'{collection_name}_no_first_gene',
            rows=sites_without_tf,
            fields=CL.RI_FIELDS
        )

        # print(rows)
        print('Finished')
    if table_name == EC.ALL_TABLES or table_name == EC.CONF_EFF_TABLE:
        collection_name = EC.TF_COLLECTION
        print(f'Generating {collection_name} collection')
        collection_ids = ids_db['identifiers']
        collection = db[collection_name]
        pipeline = [
            {
                '$project': {
                    '_id': 1,
                    'activeConformations': 1,
                    'products_ids': 1
                }
            }
        ]
        mg_objects = collection.aggregate(pipeline)
        conformations = []
        for obj in mg_objects:
            mg_tf_id = obj.get("_id", None)
            # print(mg_tf_id)
            mg_conformations = obj.get("activeConformations", None)
            if mg_conformations is None:
                conformation = {
                    '_id': mg_tf_id,
                    'type': EC.TF_COLLECTION
                }
                conformations.append(conformation)
            else:
                for conformation in mg_conformations:
                    conformations.append(conformation)
            mg_products = obj.get('products_ids')
            for product in mg_products:
                conformation = {
                    '_id': product,
                    'type': EC.PD_COLLECTION
                }
                conformations.append(conformation)
        print(f'{conformations} {len(conformations)}')


if __name__ == "__main__":
    print('Starting')
    args = arguments.load_arguments()
    id_database = args.id_database  # 'regulondbidentifiers'
    mg_database = args.database  # 'regulondbmultigenomic'
    url = args.url  # 'mongodb://localhost'
    table_name = args.table_name  # 'regulatoryInteractions'
    run(
        id_database=id_database,
        mg_database=mg_database,
        url=url,
        table_name=table_name
    )
