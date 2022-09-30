import os
import pickle
from cnvannot.common.paths import Common


def serialization_is_serialized(file_name):
    serialized_path = os.path.join(Common.serialized_path, file_name)

    if os.path.isfile(serialized_path):
        return True

    return False


def serialization_deserialize(file_name):
    with open(os.path.join(Common.serialized_path, file_name), 'rb') as f:
        return pickle.load(f)


def serialization_serialize(obj, file_name):
    with open(os.path.join(Common.serialized_path, file_name), 'wb') as f:
        pickle.dump(obj, f)
