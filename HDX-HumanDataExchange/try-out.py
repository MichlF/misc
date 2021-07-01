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
Configuration.create(hdx_site='stage', user_agent='A_Quick_Example', hdx_base_config_yaml="C:\\Users\Michl\\Documents\\GitHub\\misc\\HDX-HumanDataExchange\\.hdx_configuration.yml")
dataset = Dataset.read_from_hdx('reliefweb-crisis-app-data')
print(dataset)
print(dataset.get_dataset_date())

# dataset.set_dataset_date('2015-07-26', date_format='%Y-%m-%d')
# print(dataset.get_dataset_date())
# dataset.update_in_hdx()