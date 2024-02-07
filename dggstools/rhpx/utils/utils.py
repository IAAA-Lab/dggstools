def almost_equal(val1: float, val2: float, eps: float = 1e-6):
    return (val2 - eps) <= val1 <= (val2 + eps)


def num(s: str):
    try:
        return int(s)
    except ValueError:
        return float(s)


def insert_suffix(input_file_path: str, suffix: str) -> str:
    """
    Example: insert_suffix("a.txt","_suf") returns "a_suf.txt"
    """
    latest_dot = input_file_path.rfind(".")
    name = input_file_path[:latest_dot]
    ext_with_dot = input_file_path[latest_dot:]
    return name + suffix + ext_with_dot


def change_extension(input_file_path:str, new_ext: str) -> str:
    """
    Example: change_extension("a.txt","doc") returns "a.doc"
    """
    latest_dot = input_file_path.rfind(".")
    name = input_file_path[:latest_dot]
    return name + "." + new_ext
