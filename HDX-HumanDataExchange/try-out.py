from hdx.utilities.easy_logging import setup_logging
from hdx.hdx_configuration import Configuration
from hdx.data.dataset import Dataset

'''
See https://data.humdata.org/documentation and https://github.com/OCHA-DAP/hdx-python-api
'''

setup_logging()
# for reading data
#Configuration.create(hdx_site='prod', user_agent='A_Quick_Example', hdx_read_only=True)
# for writing data
#Configuration.create(hdx_site='stage', user_agent='A_Quick_Example', hdx_base_config_yaml="C:\\Users\Michl\\Documents\\GitHub\\private_keys\\HDX_site\\.hdx_configuration.yml")

#dataset = Dataset.read_from_hdx('reliefweb-crisis-app-data')
#print(dataset)
#print(dataset.get_dataset_date())

# dataset.set_dataset_date('2015-07-26', date_format='%Y-%m-%d')
# print(dataset.get_dataset_date())
# dataset.update_in_hdx()
#Configuration.create(hdx_site='stage', user_agent='A_Quick_Example')
Configuration.create(hdx_site='prod', user_agent='A_Quick_Example', hdx_read_only=True)
dataset = Dataset.read_from_hdx('novel-coronavirus-2019-ncov-cases')
print(dataset.get_date_of_dataset)
print(dataset.get_dataset_date())
datasets = Dataset.search_in_hdx('thailand subnational boundaries', rows=10)
print(datasets)
resources = Dataset.get_all_resources(datasets)
print(resources)
