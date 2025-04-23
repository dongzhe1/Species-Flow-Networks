import csv
import itertools
from collections import Counter
from datetime import datetime
import math
from operator import mul
from numpy import prod
from collections import defaultdict


sr = 's8'
scenarios = {
    "s1": { "compliance_rate": 0.00, "compliance_efficacy": 0.99, "noncompliance_rate": 1.00, "noncompliance_efficacy": 0.76, "nonuse_rate": 0.00, "nonuse_efficacy": 0.00, "efficacy_weighted": 0.7600},
    "s2": { "compliance_rate": 0.25, "compliance_efficacy": 0.99, "noncompliance_rate": 0.75, "noncompliance_efficacy": 0.76, "nonuse_rate": 0.00, "nonuse_efficacy": 0.00, "efficacy_weighted": 0.8175},
    "s3": { "compliance_rate": 0.50, "compliance_efficacy": 0.99, "noncompliance_rate": 0.50, "noncompliance_efficacy": 0.76, "nonuse_rate": 0.00, "nonuse_efficacy": 0.00, "efficacy_weighted": 0.8750},
    "s4": { "compliance_rate": 0.75, "compliance_efficacy": 0.99, "noncompliance_rate": 0.25, "noncompliance_efficacy": 0.76, "nonuse_rate": 0.00, "nonuse_efficacy": 0.00, "efficacy_weighted": 0.9325},
    "s5": { "compliance_rate": 1.00, "compliance_efficacy": 0.99, "noncompliance_rate": 0.00, "noncompliance_efficacy": 0.76, "nonuse_rate": 0.00, "nonuse_efficacy": 0.00, "efficacy_weighted": 0.9900},
    "s6": { "compliance_rate": 0.00, "compliance_efficacy": 0.99, "noncompliance_rate": 1.00, "noncompliance_efficacy": 0.00, "nonuse_rate": 0.00, "nonuse_efficacy": 0.00, "efficacy_weighted": 0.0000},
    "s7": { "compliance_rate": 0.00, "compliance_efficacy": 0.99, "noncompliance_rate": 0.75, "noncompliance_efficacy": 0.76, "nonuse_rate": 0.25, "nonuse_efficacy": 0.00, "efficacy_weighted": 0.5700},
    "s8": { "compliance_rate": 0.25, "compliance_efficacy": 0.99, "noncompliance_rate": 0.25, "noncompliance_efficacy": 0.76, "nonuse_rate": 0.50, "nonuse_efficacy": 0.00, "efficacy_weighted": 0.4375},
    "baseline": {"compliance_rate": 1, "compliance_efficacy": 0, "noncompliance_rate": 0, "noncompliance_efficacy": 0.76, "nonuse_rate": 0, "nonuse_efficacy": 0.00, "efficacy_weighted": 0},
}

# %%
def GetPortData(fn, field, delim):
    ports = {}
    with open(fn) as csvfile:
        reader = csv.DictReader(csvfile, delimiter=delim)
        for row in reader:
            ports[row[field]] = row
    return ports


def filter_empty_temp_sal(ports):
    for key, val in list(ports.items()):
        if ((ports[key]['YR_MEAN_T'] == 'NA') or (ports[key]['Salinity'] == 'NA')):
            ports.pop(key)
    return ports


def FilterByEnv(ports, pair, TempThreshold, SalThreshold):
    FromPort, ToPort = pair
    if abs(float(ports[FromPort]['YR_MEAN_T']) - float(ports[ToPort]['YR_MEAN_T'])) >= TempThreshold or abs(
            float(ports[FromPort]['Salinity']) - float(ports[ToPort]['Salinity'])) >= SalThreshold:
        return False
    else:
        return True


def BuildDIPenv(TempTolerances, SalTolerances, ports, DIP):
    DIPE = {}
    for tt in TempTolerances:
        for st in SalTolerances:
            DIPE[tt, st] = set([pair for pair in DIP if FilterByEnv(ports, pair, tt, st)])
    return DIPE


def ReadShipFreqs(fn):
    edges = {}
    with open(fn) as f:
        for line in f:
            fields = line.strip().split(',')
            FromNode = fields[0]
            ToNode = fields[1]
            weight = int(fields[2])
            edges[FromNode, ToNode] = weight
    return edges


def fouling(source, dest, duration, stay_d, distance, antifouling_p):
    v = distance / duration
    if (abs(float(ports[source]['LATITUDE_DECIMAL'])) < 0.35):  # trpoical
        f_risk = (0.000000129 * stay_d ** 3 - 0.000083165 * stay_d ** 2 + 0.01495187 * stay_d) * (
            antifouling_p) * math.exp(-0.008 * v)

    else:  # temperate
        f_risk = (0.0000000014 * stay_d ** 3 - 0.000016566 * stay_d ** 2 + 0.00519377 * stay_d) * (
            antifouling_p) * math.exp(-0.008 * v)
    return f_risk


def deg2rad(deg):
    return deg * (math.pi / 180)


def getDistanceFromLatLonInKm(source, dest):
    # print([lat1,lon1,lat2,lon2])
    R = 6371;  # Radius of the earth in km

    lat1 = float(ports[source]['LATITUDE_DECIMAL'])
    lon1 = float(ports[source]['LONGITUDE_DECIMAL'])
    lat2 = float(ports[dest]['LATITUDE_DECIMAL'])
    lon2 = float(ports[dest]['LONGITUDE_DECIMAL'])

    dLat = deg2rad(lat2 - lat1)  # deg2rad below
    dLon = deg2rad(lon2 - lon1);
    a = math.sin(dLat / 2) * math.sin(dLat / 2) + math.cos(deg2rad(lat1)) * math.cos(deg2rad(lat2)) * math.sin(
        dLon / 2) * math.sin(dLon / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = R * c;  # Distance in km
    return d;


# only if they are both from the same ecoregion the chance of being alien is 0
# if source is in neighboring ecoregion of dest, chance of being alien is half.
# OW. chance of being alien is 1.
def filterbyeco(source, dest):
    p_alien = 1
    s_eco = ports[source]['MEOW_region']
    d_eco = ports[dest]['MEOW_region']
    if ((s_eco != 'NA') and (d_eco != 'NA')):
        if (s_eco == d_eco):
            p_alien = 0

        elif ((s_eco != 'NA') and (s_eco in ports[dest]['MEOW_Neighbors'].split('|'))):
            p_alien = 0
        else:
            p_alien = 1

    if ((s_eco == 'NA') and (
            d_eco != 'NA')):  # source is a freshwater port. only need to check the freshwater neighbors of
        # dest, since only in that case species will be native.
        s_eco = ports[source]['FEOW_region']
        if ((s_eco != 'NA') and ((s_eco in ports[dest]['FEOW_Neighbors'].split('|')) or (
                s_eco in ports[dest]['MEOW_Neighbors'].split('|')))):
            p_alien = 0
        else:
            p_alien = 1

    if ((d_eco == 'NA') and (s_eco != 'NA')):  # dest is a freshwater port. only need to check the fresh neighbors of
        # dest, since only in that case species will be native.
        d_eco = ports[dest]['FEOW_region']
        if ((d_eco != 'NA') and ((s_eco in ports[dest]['FEOW_Neighbors'].split('|')) or (
                s_eco in ports[dest]['MEOW_Neighbors'].split('|')))):
            p_alien = 0
        else:
            p_alien = 1

    return (p_alien)


def establishment(source, dest):
    try:
        sal_diff = abs(float(ports[source]['Salinity']) - float(ports[dest]['Salinity']))
        temp_diff = abs(float(ports[source]['YR_MEAN_T']) - float(ports[dest]['YR_MEAN_T']))
        p_estab = 0.00015 * math.exp(
            -0.5 * ((temp_diff / 2) ** 2 + (sal_diff / 10) ** 2))  # sigma_t=2'C & sigma_s=10 ppt
        return (p_estab)
    except:
        return (-1)


def IsValid_Month(pair):
    first = pair[0]
    second = pair[1]
    FirstSailYear = first['SAIL DATE'].split(' ')[0][0:4]
    SecondArrivalYear = second['ARRIVAL DATE'].split(' ')[0][0:4]
    # same vessel
    if first['VESSEL ID'] == second['VESSEL ID'] and FirstSailYear != '' and SecondArrivalYear != '':
        if not first['PLACE ID'] == second['PLACE ID']:  # different ports
            if FirstSailYear == SecondArrivalYear or FirstSailYear in (
            '1997', '1999', '2002', '2005', '2008', '2012', '2015', '2018') and (
                    int(SecondArrivalYear) == int(FirstSailYear) + 1):
                try:
                    trip_durtaion = float(second['DURATION'])
                    stay_at_source_duration = float(second['STAY DURATION'])
                    if (trip_durtaion > 0 and stay_at_source_duration >= 0):
                        return True, trip_durtaion, stay_at_source_duration
                except:
                    print('error')
                    return False, -1, -1
    return False, -1, -1


def filterbyeco_same_only(source, dest):
    try:
        p_alien = 1
        s_eco = ports[source]['MEOW_region']
        d_eco = ports[dest]['MEOW_region']

        if ((s_eco != 'NA') and (d_eco != 'NA')):
            if (s_eco == d_eco):
                p_alien = 0
            else:
                p_alien = 1

        return (p_alien)
    except:
        return (-1)


def DumpEdges(edges, fn):
    with open(fn, 'w') as f:
        for edge in edges:
            FromPort, ToPort = edge
            f.write(','.join(map(str,
                                 [FromPort, ToPort,
                                  ports[FromPort]['NAME'], ports[ToPort]['NAME'],
                                  ports[FromPort]['LATDD'], ports[FromPort]['LONGDD'],
                                  ports[ToPort]['LATDD'], ports[ToPort]['LONGDD'],
                                  edges[edge],
                                  ports[FromPort]['YR_MEAN_T'], ports[ToPort]['YR_MEAN_T'],
                                  ports[FromPort]['Salinity'], ports[ToPort]['Salinity']]))
                    + '\n')


def DumpWeightedEdges(edges, network, fn):
    with open(fn, 'w') as f:
        for edge in edges:
            try:
                FromPort, ToPort = edge
                f.write(','.join(map(str,
                                     [FromPort, ToPort,
                                      network[edge],
                                      ports[FromPort]['NAME'], ports[ToPort]['NAME'],
                                      ports[FromPort]['ABBREV'], ports[ToPort]['ABBREV'],
                                      ports[FromPort]['LATDD'], ports[FromPort]['LONGDD'],
                                      ports[ToPort]['LATDD'], ports[ToPort]['LONGDD'],
                                      edges[edge],
                                      ports[FromPort]['YR_MEAN_T'], ports[ToPort]['YR_MEAN_T'],
                                      ports[FromPort]['Salinity'], ports[ToPort]['Salinity']]))
                        + '\n')
            except Exception as e:
                print(e)


def DumpPorts(d, fn):
    with open(fn, 'w') as f:
        for port in d:
            try:
                f.write(','.join(map(str,
                                     [port, d[port],
                                      ports[port]['NAME'], ports[port]['ABBREV'],
                                      ports[port]['LATDD'], ports[port]['LONGDD'],
                                      ports[port]['YR_MEAN_T'], ports[port]['Salinity'],
                                      ]))
                        + '\n')
            except Exception as e:
                print(e)

def calculate_p_intro(move_discharge, trip_duration):
    ballast_discharge = (1 - math.exp(-3.22 * (10**-6) * move_discharge))
    scenario = scenarios[sr]
    efficacy_ratio = (scenario["compliance_rate"] * (1 - scenario["compliance_efficacy"])
                      + scenario["noncompliance_rate"] * (1 - scenario["noncompliance_efficacy"])
                      + scenario["nonuse_rate"] * (1 - scenario["nonuse_efficacy"]))
    p_intro = efficacy_ratio * ballast_discharge * math.exp(-0.02 * trip_duration)
    return p_intro

year = '2018'
y = year + '/'
order = 2

env='noEnv_'
eco = 'noEco_'
paul = ''
clean_move = '../../data/moves/moves_cleaned_' + str(year) + '.txt'
port_data = '../../data/Places_allportdata_mergedSept2017.csv'

# output
InputForHON_Foulig = '../../data/' + y + paul + 'Fouling_' + env + eco + sr + '.csv'
InputForHON_Ballast = '../../data/' + y + paul + 'Ballast_' + env + eco + sr + '.csv'

# %%
ports = GetPortData(port_data, 'ID', ',')
print("Number of all ports: ", len(ports))

# %%
moves = []
with open(clean_move) as csvfile:
    reader = csv.DictReader(csvfile, delimiter='|')
    for row in reader:
        moves.append(row)
print("Number of all moves: ", len(moves))


ShipID = -1
ShipTraj = []
ProbsRaw = {}
ErrorCounter_0 = 0
for prevline in range(len(moves) - 1):

    nextline = prevline + 1
    prevmove, nextmove = moves[prevline], moves[nextline]
    pair = (prevmove, nextmove)
    source, dest = prevmove['PLACE ID'], nextmove['PLACE ID']
    PortPair = (source, dest)
    IsValid, trip_duration, stay_duration = IsValid_Month(pair)

    if IsValid:
        if prevmove['VESSEL ID'] != ShipID:  # if different ships than the last:
            ShipTraj = [source, dest]
        else:
            ShipTraj.append(nextmove['PLACE ID'])

        ShipID = prevmove['VESSEL ID']
        try:
            distance = getDistanceFromLatLonInKm(source, dest)
            subseq = []  # subsequences
            for lastn in range(2, 1 + min(len(ShipTraj), order)):
                subseq.append(ShipTraj[-lastn:])

            if not PortPair in ProbsRaw:
                ProbsRaw[PortPair] = []
            ProbsRaw[PortPair].append({'ship': prevmove['VESSEL ID'],
                                       'trip_duration': float(trip_duration),
                                       'stay_duration': float(stay_duration),
                                       'ballast_discharge': float(prevmove['BALLAST DISCHARGE']),
                                       'antifouling_p': float(prevmove['ANTIFOULING PROB']),
                                       'GWT': prevmove['GWT'],
                                       'distance': float(distance),
                                       'subseq': subseq})
        # else:
        except:
            ErrorCounter_0 += 1
    else:
        # reset cached ship trajectory
        # if a ship sails out of the Arctic and comes back in, treat as two trajectories.
        ShipID = -1
print("ProbsRaw Made")
print(ErrorCounter_0)
len(ProbsRaw)

# %%
####################################################################
# Ballast & Fouling
####################################################################
ballast_w, foul_w = [], []
counter = 0;
ErrorCounter = 0
for pair in ProbsRaw:
    for move in ProbsRaw[pair]:
        if eco == 'noEco_':
            p_alien = 1
        elif eco == 'sameEco_':
            p_alien = filterbyeco_same_only(pair[0], pair[1])
        else:
            p_alien = filterbyeco(pair[0], pair[1])

        if env == 'noEnv_':
            p_estab = 1
        else:
            p_estab = establishment(pair[0], pair[1])

        p_intro = calculate_p_intro(move['ballast_discharge'], trip_duration)
        if (p_estab != -1):  # portss who doesn't have env data are ignored.
            prob_ballast = p_alien * p_intro * p_estab
            ballast_w.append(prob_ballast)

            p_fouling = fouling(pair[0], pair[1], move['trip_duration'], move['stay_duration'], move['distance'],
                                move['antifouling_p'])
            prob_fouling = p_alien * p_fouling * p_estab
            foul_w.append(prob_fouling)

print(counter, " errors")
max_ballast = max(ballast_w)
max_fouling = max(foul_w)
print('max_ballast: ', max_ballast)
print('max_fouling: ', max_fouling)

probs_ballast = defaultdict(list)
probs_fouling = defaultdict(list)
b = open(InputForHON_Ballast, 'w')
#f = open(InputForHON_Foulig, 'w')

for pair in ProbsRaw:
    for move in ProbsRaw[pair]:
        if (True):
            if eco == 'noEco_':
                p_alien = 1
            else:
                p_alien = filterbyeco_same_only(pair[0], pair[1])
            if env == 'noEnv_':
                p_estab = 1
            else:
                p_estab = establishment(pair[0], pair[1])

            p_intro = calculate_p_intro(move['ballast_discharge'], trip_duration)

            if (p_estab != -1):  # portss who doesn't have env data are ignored.
                prob_ballast = p_alien * p_intro * p_estab
                ballast_w.append(prob_ballast)
                # if(prob_ballast>0):
                probs_ballast[pair].append(prob_ballast)
                for subseq in move['subseq']:
                    b.write(' '.join(subseq) + ' ' + str(prob_ballast) + '\n')
                    break

                if False:
                    p_fouling = fouling(pair[0], pair[1], move['trip_duration'], move['stay_duration'], move['distance'],
                                        move['antifouling_p'])
                    prob_fouling = p_alien * p_fouling * p_estab / max_fouling
                    foul_w.append(prob_fouling)
                    probs_fouling[pair].append(prob_fouling)
                    for subseq in move['subseq']:
                        #f.write(' '.join(subseq) + ' ' + str(prob_fouling) + '\n')
                        break

        else:
            ErrorCounter += 1

print(ErrorCounter, " biofouling errors")
#f.close()
b.close()
print('Written InputForHoN ')
