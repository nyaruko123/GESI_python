"""
Table of correspondence between HL and SPL
Function [Table] = SPLatHL0dB_Table()
    OUTPUT:
        Table['freq']: Frequency
        Table['SPLatHL0dB']: SPL in dB
"""


def SPLatHL0dB_Table():
    # Table including intermediate frequencies from ANSI_S3.6-2010 (1996)
    Freq_SPLdBatHL0dB_List = [
        (125, 45.0), (160, 38.5), (200, 32.5),
        (250, 27.0), (315, 22.0), (400, 17.0),
        (500, 13.5), (630, 10.5), (750, 9.0), (800, 8.5),
        (1000, 7.5), (1250, 7.5), (1500, 7.5), (1600, 8.0),
        (2000, 9.0), (2500, 10.5), (3000, 11.5), (3150, 11.5),
        (4000, 12.0), (5000, 11.0), (6000, 16.0), (6300, 21.0),
        (8000, 15.5)
    ]

    Table = {
        'freq': [pair[0] for pair in Freq_SPLdBatHL0dB_List],
        'SPLatHL0dB': [pair[1] for pair in Freq_SPLdBatHL0dB_List],
        'Speech': 20.0,
        'Standard': 'ANSI-S3.6_2010',
        'Earphone': 'Any supra aural earphone having the characteristics described in clause 9.1.1 or ISO 389-1',
        'ArtificialEar': 'IEC 60318-1'
    }

    return Table


