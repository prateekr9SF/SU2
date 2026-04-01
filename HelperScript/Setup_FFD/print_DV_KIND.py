"""
Generates a string that repeats "FFD_CONTROL_POINT" 10 times,
separated by commas, and appends "FFD_ROTATION" at the end.
"""

def generate_ffd_string(count):
    return ", ".join(["FFD_CONTROL_POINT"] * count) + ", FFD_ROTATION"

def print_decimal_values(count):
    print(", ".join(["0.001"] * count))

if __name__ == "__main__":
    

    count = 206


    
    print(generate_ffd_string(count))
    
    print("\n")
    
    print_decimal_values(count+1)
    

