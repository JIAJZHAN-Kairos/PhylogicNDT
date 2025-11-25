# from enum import Enum, unique #requires >2.7.6 or 3.x

def Enum(**enums):
    return type('Enum', (), enums)


Cluster = Enum(**dict([('C1', 'green'), ('C2', 'blue'), ('C3', 'red')]))

# csize contains chromosome bp lengths
CSIZE = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
            133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
            58617616, 64444167, 46709983, 50818468, 156040895, 57227415
        ]

chrom_dict = {'chr' + str(c): s for c, s in zip(list(range(1, 23)) + ['X', 'Y'], CSIZE)}
ChromSize = Enum(**chrom_dict)

# centromeres (define arm-level lengths)
CENT_LOOKUP = {1: 123400000, 2: 93900000, 3: 90900000, 4: 50000000, 5: 48800000,
               6: 59800000, 7: 60100000, 8: 45200000, 9: 43000000, 10: 39800000,
               11: 53400000, 12: 35500000, 13: 17700000, 14: 17200000, 15: 19000000,
               16: 36800000, 17: 25100000, 18: 18500000, 19: 26200000, 20: 28100000, 21: 12000000,
               22: 15000000, 23: 61000000, 24: 10400000}

MutStatus = Enum(OK="OK",
                 REMOVED="REMOVED",  # blacklisted
                 GRAYLIST="GRAYLIST")  # not used in clustering

MutType = Enum(INS="INS",
               DEL="DEL",
               SNV="SNV")
