import numpy as np
import pyparsing as pa
import matplotlib.pyplot as plt

# enable catching 
#pa.ParserElement.enablePackrat()

# Parser
floatNumber = pa.Regex(r'(\-)?\d+(\.)(\d*)?([eE][\-\+]\d+)?')
natural = pa.Word(pa.nums)

def parse_file(file_name):
    parser = generate_parser()
    xs = parser.parseFile(file_name).asList()

    return np.array(xs, dtype=np.float)

def parse_section(label, number=2, number_parser=floatNumber):
    return pa.Suppress(pa.SkipTo(label)) + pa.Suppress(pa.SkipTo(number_parser)) + number_parser * number

def generate_parser():
    """
    This parser takes structural info from the output of a cp2k calculation
    """
    parser_step = parse_section('STEP NUMBER', number=1, number_parser=natural)
    parser_pressure = parse_section(' PRESSURE [bar]')
    parser_volumen =  parse_section(' VOLUME[bohr^3]')
    parser_cell_lnth = parse_section('AVE. CELL LNTHS[bohr]', 3)
    parser_cell_angl = parse_section( 'AVE. CELL ANGLS[deg]', 3)

    parser_section = parser_step + parser_pressure + parser_volumen + parser_cell_lnth + parser_cell_angl

    return pa.OneOrMore(pa.Group(parser_section))



if __name__ == "__main__":
    parse_file('test.out')


