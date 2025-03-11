"""
 Title:         Helper
 Description:   General helper functions
 References:    https://stackoverflow.com/questions/8391411/how-to-block-calls-to-print
 Author:        Janzen Choi

"""

# Libraries
import math, os, sys
import numpy as np

def periodify(value:float, min:float, max:float) -> float:
    """
    Applies periodicity to a value

    Parameters:
    * `value`: The value to be applied
    * `min`:   The lower bound of the periodicity
    * `max`:   The upper bound of the periodicity
    
    Returns the periodified value
    """
    if value >= min and value <= max:
        return value
    range = max - min
    if value < min:
        return periodify(value+range, min, max)
    if value > max:
        return periodify(value-range, min, max)

def remove_nan(data_list:list) -> list:
    """
    Removes nan values from a list of data values

    Parameters:
    * `data_list`: The list of data values

    Returns the list of data values without nan values
    """
    return list(filter(lambda x: not math.isnan(x), data_list))

def get_closest(x_list:list, y_list:list, x_value:float) -> float:
    """
    Finds the closest corresponding y value given an x value;
    does not interpolate

    Parameters:
    * `x_list`:  The list of x values
    * `y_list`:  The list of y values
    * `x_value`: The x value to get the closest value of
    
    Returns the closest value
    """
    x_diff_list = [abs(x-x_value) for x in x_list]
    x_min_diff = min(x_diff_list)
    x_min_index = x_diff_list.index(x_min_diff)
    return y_list[x_min_index]

def quick_spline(x_list:list, y_list:list, x_value:float) -> float:
    """
    Conducts a quick evaluation using spline interpolation without
    conducting the whole interpolation; assumes that the x_value is
    between min(x_list) and max(x_list) and that x_list is sorted

    Parameters:
    * `x_list`:  The list of x values
    * `y_list`:  The list of y values
    * `x_value`: The x value to evaluate
    
    Returns the evaluated y value
    """
    if len(x_list) != len(y_list):
        raise ValueError("Length of lists do not match!")
    for i in range(len(x_list)-1):
        if x_list[i] <= x_value and x_value <= x_list[i+1]:
            gradient = (y_list[i+1]-y_list[i])/(x_list[i+1]-x_list[i])
            y_value = gradient*(x_value - x_list[i]) + y_list[i]
            return y_value
    return None

def sort_dict(unsorted_dict:dict) -> dict:
    """
    Sorts a dictionary in ascending order based on its values
    
    Parameters:
    * `unsorted_dict`: Dictionary with unsorted values
    
    Returns the sorted dictionary
    """
    sorted_dict = dict(sorted(unsorted_dict.items(), key=lambda item: item[1]))
    return sorted_dict

def remove_consecutive_duplicates(value_list:list) -> list:
    """
    Removes duplicates from a list of values with the
    assumption that the duplicates are consecutive;
    maintains order
    """
    new_value_list = []
    for value in value_list:
        if not value in new_value_list:
            new_value_list.append(value)
    return new_value_list

def flatten(list_of_lists:list) -> list:
    """
    Flattens a 2D list into a 1D list
    
    Parameters:
    * `list_of_lists`: A list of lists (i.e., a 2D grid)
    
    Returns the flattened list
    """
    return [item for sublist in list_of_lists for item in sublist]

def try_float(value:str) -> float:
    """
    Tries to float cast a value

    Parameters:
    * `value`: The value to be float casted

    Returns the value
    """
    try:
        value = float(value)
    except:
        pass
    return value

def transpose(list_of_lists:list) -> list:
    """
    Transposes a 2D list of lists
    
    Parameters:
    * `list_of_lists`: A list of lists (i.e., a 2D grid)
    
    Returns the transposed list of lists
    """
    transposed = np.array(list_of_lists).T.tolist()
    return transposed

def integer_to_ordinal(n:int):
    """
    Converts an integer to an ordinal string
    """
    if 10 <= n % 100 <= 20:
        suffix = "th"
    else:
        suffix = {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return str(n) + suffix

def get_thinned_list(unthinned_list:list, density:int) -> list:
    """
    Gets a thinned list

    Parameters:
    * `unthinned_list`: The list before thinning
    * `density`:        The goal density of the thinned list

    Returns the thinned list
    """
    src_data_size = len(unthinned_list)
    step_size = src_data_size / density
    thin_indexes = [math.floor(step_size*i) for i in range(1, density - 1)]
    thin_indexes = [0] + thin_indexes + [src_data_size - 1]
    thinned_list = [unthinned_list[i] for i in thin_indexes]
    return thinned_list

def get_sorted(value_list:list, reverse:bool=True) -> tuple:
    """
    Gets the top values and indexes of a list of values
    
    Parameters:
    * `value_list`: The list of values
    
    Returns the list of top values and indexes
    """
    sorted_value_list = sorted(value_list, reverse=reverse)
    sorted_index_list = []
    for value in sorted_value_list:
        for i in range(len(value_list)):
            if value == value_list[i] and not i in sorted_index_list:
                sorted_index_list.append(i)
                break
    return sorted_value_list, sorted_index_list

def pad_to_length(value_str:str, length:int) -> str:
    """
    Pads a stringified number with zeroes to match the desired length
    
    Parameters:
    * `value_str`: The stringified number
    * `length`:    The desired length

    Returns the padded stringified number
    """
    value_str = str(value_str) # force string
    curr_length = len(value_str)
    if curr_length < length:
        difference = length - curr_length
        return value_str + "0"*difference
    return value_str # no changes

def round_sf(value:float, sf:int) -> float:
    """
    Rounds a float to a number of significant figures

    Parameters:
    * `value`: The value to be rounded; accounts for lists
    * `sf`:    The number of significant figures

    Returns the rounded number
    """
    if isinstance(value, list):
        return [round_sf(v, sf) for v in value]
    format_str = "{:." + str(sf) + "g}"
    rounded_value = float(format_str.format(value))
    return rounded_value

class BlockPrint:
    """
    Blocks print messages
    """

    def __enter__(self) -> None:
        """
        Auxiliary function
        """
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """
        Auxiliary function
        """
        sys.stdout.close()
        sys.stdout = self._original_stdout
