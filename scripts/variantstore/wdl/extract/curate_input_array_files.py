import numpy as np
from contextlib import contextmanager

SAMPLE_MAP_FILE_SUFFIX = "sample_map_file"
SAMPLE_NAME_FILE_SUFFIX = "sample_name_list_file"
VCF_FILE_SUFFIX = "vcf_list_file"
VCF_INDEX_FILE_SUFFIX = "vcf_index_list_file"

@contextmanager
def handle_file_error(file_name):
    try:
        yield
    except:
        print(f"ERROR: required file named '{file_name}' does not exist.")


def curate_input_arrays():
    sample_map_array = vcf_array = vcf_indexes_array = sample_names_array = []
    with handle_file_error(f"input_{SAMPLE_MAP_FILE_SUFFIX}"):
        sample_map_array = np.loadtxt(f"input_{SAMPLE_MAP_FILE_SUFFIX}", dtype=str, delimiter=",")
    with handle_file_error(f"input_{VCF_FILE_SUFFIX}"):
        vcf_array = np.loadtxt(f"input_{VCF_FILE_SUFFIX}", dtype=str)
    with handle_file_error(f"input_{VCF_INDEX_FILE_SUFFIX}"):
        vcf_indexes_array = np.loadtxt(f"input_{VCF_INDEX_FILE_SUFFIX}", dtype=str)
    with handle_file_error(f"input_{SAMPLE_NAME_FILE_SUFFIX}"):
        sample_names_array = np.loadtxt(f"input_{SAMPLE_NAME_FILE_SUFFIX}", dtype=str)
    rows_to_delete = []

    # use input_sample_names_array to figure out which index "rows" to delete
    for i in range(len(sample_names_array)):
        if sample_names_array[i] not in sample_map_array:
            rows_to_delete.append(i)

    # re-create input arrays using array of "rows" to delete
    vcf_array = [vcf_array[i] for i in range(len(vcf_array)) if i not in rows_to_delete]
    vcf_indexes_array = [vcf_indexes_array[i] for i in range(len(vcf_indexes_array)) if
                         i not in rows_to_delete]
    sample_names_array = [sample_names_array[i] for i in range(len(sample_names_array)) if
                          i not in rows_to_delete]

    print(f"Creating 'output_{SAMPLE_NAME_FILE_SUFFIX}', 'output_{VCF_FILE_SUFFIX}' and 'output_{VCF_INDEX_FILE_SUFFIX}'.")
    np.savetxt(f"output_{SAMPLE_NAME_FILE_SUFFIX}", sample_names_array, fmt='%s')
    np.savetxt(f"output_{VCF_FILE_SUFFIX}", vcf_array, fmt='%s')
    np.savetxt(f"output_{VCF_INDEX_FILE_SUFFIX}", vcf_indexes_array, fmt='%s')


if __name__ == '__main__':
    curate_input_arrays()
