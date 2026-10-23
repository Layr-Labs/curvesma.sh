#!/usr/bin/env python3
"""Exact-rational empirical tables; floating fields are explanatory displays only.
No statistic in this file is a proved production mixing or probability bound.
"""
from collections import Counter
from fractions import Fraction
import json
import math
from pathlib import Path

GRID = [(1,2),(1,1),(3,2),(2,1),(3,1)]


def rational(x):
    x = Fraction(x)
    return {'numerator':str(x.numerator), 'denominator':str(x.denominator)}


def grid_count(items, N, a, b, key):
    return sum(row[key]**2*b*b <= a*a*N for row in items)


def summarize(data):
    if 'pipeline_file' in data:
        assert data['pipeline_file'] == 'mechanism-pipeline.json'
        data = {**data, 'pipeline': json.loads(Path(__file__).with_name('mechanism-pipeline.json').read_text())['pipeline']}
    total = sum(len(f['first_cycles']) for f in data['fields'])
    cdf, fields, pipes = [], [], []
    for key in ('rho','lambda','mu'):
        for a,b in GRID:
            actual = sum(grid_count(f['first_cycles'],f['N'],a,b,key) for f in data['fields'])
            ref = sum(grid_count(f['reference_cycles'],f['N'],a,b,key) for f in data['fields'])
            cdf.append({'quantity':key, 'threshold':f'{a}/{b}', 'n':total,
                        'concrete':actual, 'random_function':ref,
                        'absolute_difference':rational(Fraction(abs(actual-ref),total))})
    for fi, f in enumerate(data['fields']):
        n, N = len(f['first_cycles']), f['N']
        # Finite table occupancy statistic R*sum(c_j^2)/m - m is exact.
        m = f['uniform_coordinate_samples']
        chi = Fraction(1024*sum(v*v for v in f['partition_counts']),m)-m
        lambdas = sorted(t['lambda'] for t in f['first_cycles'])
        fields.append({'index':fi,'p':f['p'],'r':f['r'],'N':N,
                       'log2_N_display':math.log2(N),'independent_tables':n,
                       'min_lambda':min(lambdas),'median_lambda':lambdas[n//2],
                       'max_lambda':max(lambdas),
                       'min_rho':min(t['rho'] for t in f['first_cycles']),
                       'max_rho':max(t['rho'] for t in f['first_cycles']),
                       'cycles_lambda_le_32':sum(v<=32 for v in lambdas),
                       'cycles_lambda_le_1024':sum(v<=1024 for v in lambdas),
                       'partition_pearson_chi_square':rational(chi),'partition_df':1023,
                       'uniform_samples':m,'uniform_mark_counts':f['uniform_mark_counts'],
                       'uniform_density_ratio': [rational(Fraction(v*(1<<d),m))
                                                 for v,d in zip(f['uniform_mark_counts'],(8,10,12))],
                       'first_path_distinct_states':sum(t['rho'] for t in f['first_cycles']),
                       'first_path_mark_counts':[sum(t['marks'][i] for t in f['first_cycles']) for i in range(3)],
                       'first_path_mark_free_cycles':[sum(t['cycle_marks'][i]==0 for t in f['first_cycles']) for i in range(3)],
                       'long_basin_trials':len(f['basin_trials']),
                       'cycles_per_long_trial':[len(t['cycles']) for t in f['basin_trials']],
                       'basin_mark_free_cycles':[sum(c['cycle_marks'][i]==0 for t in f['basin_trials'] for c in t['cycles']) for i in range(3)],
                       'lambda_histogram':dict(sorted(Counter(lambdas).items()))})
    for fi in range(len(data['fields'])):
        for d in (8,10,12):
            group = [p for p in data['pipeline'] if p['field_index']==fi and p['density_bits']==d]
            if not group:
                continue
            agg = {k:sum(p[k] for p in group) for k in ('windows','tortoise_iterations','evaluations',
                    'dense_multiplication_analogues','starts','reports','floyd_detections','timeout_restarts',
                    'cycle_basin_iterations','verified_floyd_relations','verified_report_relations',
                    'degenerate_floyd_relations','successful_windows')}
            T, N = agg['tortoise_iterations'], group[0]['N']
            collisions = agg['verified_floyd_relations']+agg['verified_report_relations']
            agg.update(field_index=fi, density_bits=d,
                       success_fraction=rational(Fraction(agg['successful_windows'],agg['windows'])),
                       cycle_basin_loss_fraction=rational(Fraction(agg['cycle_basin_iterations'],T)),
                       observed_mark_density_ratio=rational(Fraction(agg['reports']*(1<<d),T)),
                       verified_collisions_per_sqrt_R_normalized_budget=rational(Fraction(32*collisions,T)),
                       verified_collisions_per_sqrt_N_normalized_budget_display=collisions*math.sqrt(N)/T,
                       budget_per_window_over_sqrt_N_display=group[0]['budget_per_window']/math.sqrt(N),
                       table_success_fractions=[rational(Fraction(p['successful_windows'],p['windows'])) for p in group])
            pipes.append(agg)
    output = {'seed':data['seed'],'elapsed_seconds':data['elapsed_seconds'],
              'cycles_observed':data['cycles_observed'],'independent_table_first_cycles':total,
              'fields':fields,'empirical_CDF':cdf,'pipeline':pipes,
              'empirical_grid_deviation': {key:rational(max(Fraction(int(row['absolute_difference']['numerator']),int(row['absolute_difference']['denominator'])) for row in cdf if row['quantity']==key)) for key in ('rho','lambda','mu')},
              'interpretation':'Grid deviations compare empirical CDFs over independent tables. Restart samples within one table are basin-weighted diagnostics, not independent map samples. Neither these deviations nor two-cycle exclusion is a production statistical-distance bound. R=1024 is the jump-table width; sqrt_R normalization is separate from the birthday scale sqrt_N. Pipeline counts are detected/recovered relations, not independent birthday events.'}
    return output


if __name__ == '__main__':
    data = json.loads(Path(__file__).with_name('mechanism-output.json').read_text())
    print(json.dumps(summarize(data),indent=2))
