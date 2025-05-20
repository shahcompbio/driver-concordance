import pytest

import numpy as np
import pandas as pd

from gdan_hcmi_tools.process import extract_sample_ids, interp_nan, interp_nan_per_chrom, Cohort, transform_aa_change

@pytest.fixture
def links():
    return pd.DataFrame([
        ['C1', 'M1', 'T1', 'N1', '',],
        ['C2', 'M2', 'T2', 'N2', 'ME2',],
    ], columns=['case', 'model', 'tumor', 'normal', 'model_expanded'])

def test_cohort_1(links): # init test
    cohort = Cohort(links)
    assert set(cohort.tumor_ids) == {'T1', 'T2'}, cohort.tumor_ids
    assert set(cohort.model_ids) == {'M1', 'M2', 'ME2'}, cohort.model_ids
    assert set(cohort.sample_ids) == {'T1', 'T2', 'M1', 'M2', 'ME2'}, cohort.sample_ids
    assert cohort.model2tumor == {
        'M1': 'T1',
        'M2': 'T2',
        'ME2': 'T2',
    }, cohort.model2tumor
    assert set(cohort.tumor2models.keys()) == {'T1', 'T2'}, cohort.tumor2models
    assert set(cohort.tumor2models['T1']) == {'M1',}, cohort.tumor2models
    assert set(cohort.tumor2models['T2']) == {'M2', 'ME2',}, cohort.tumor2models

def test_cohort_2(links): # exclude T1
    cohort = Cohort(links)
    exclude = ['T1']
    cohort.tumor_ids = [t for t in cohort.tumor_ids if t not in exclude]
    assert set(cohort.tumor_ids) == {'T2'}, cohort.tumor_ids
    assert set(cohort.model_ids) == {'M2', 'ME2'}, cohort.model_ids
    assert set(cohort.sample_ids) == {'T2', 'M2', 'ME2'}, cohort.sample_ids
    assert cohort.model2tumor == {
        'M2': 'T2',
        'ME2': 'T2',
    }, cohort.model2tumor
    assert set(cohort.tumor2models.keys()) == {'T2'}, cohort.tumor2models
    assert set(cohort.tumor2models['T2']) == {'M2', 'ME2',}, cohort.tumor2models

def test_cohort_3(links): # exclude T2
    cohort = Cohort(links)
    exclude = ['T2']
    cohort.tumor_ids = [t for t in cohort.tumor_ids if t not in exclude]
    assert set(cohort.tumor_ids) == {'T1'}, cohort.tumor_ids
    assert set(cohort.model_ids) == {'M1'}, cohort.model_ids
    assert set(cohort.sample_ids) == {'T1', 'M1'}, cohort.sample_ids
    assert cohort.model2tumor == {
        'M1': 'T1',
    }, cohort.model2tumor
    assert set(cohort.tumor2models.keys()) == {'T1'}, cohort.tumor2models
    assert set(cohort.tumor2models['T1']) == {'M1'}, cohort.tumor2models

def test_cohort_4(links): # regenerate original after T1 removal
    cohort = Cohort(links)
    exclude = ['T1']
    cohort.tumor_ids = [t for t in cohort.tumor_ids if t not in exclude]
    cohort.tumor_ids = ['T1', 'T2']    
    assert set(cohort.tumor_ids) == {'T1', 'T2'}, cohort.tumor_ids
    assert set(cohort.model_ids) == {'M1', 'M2', 'ME2'}, cohort.model_ids
    assert set(cohort.sample_ids) == {'T1', 'T2', 'M1', 'M2', 'ME2'}, cohort.sample_ids
    assert cohort.model2tumor == {
        'M1': 'T1',
        'M2': 'T2',
        'ME2': 'T2',
    }, cohort.model2tumor
    assert set(cohort.tumor2models.keys()) == {'T1', 'T2'}, cohort.tumor2models
    assert set(cohort.tumor2models['T1']) == {'M1',}, cohort.tumor2models
    assert set(cohort.tumor2models['T2']) == {'M2', 'ME2',}, cohort.tumor2models

def test_cohort_5(links): # exclude M1
    cohort = Cohort(links)
    exclude = ['M1']
    cohort.model_ids = [t for t in cohort.model_ids if t not in exclude]
    assert set(cohort.tumor_ids) == {'T2'}, cohort.tumor_ids
    assert set(cohort.model_ids) == {'M2', 'ME2'}, cohort.model_ids
    assert set(cohort.sample_ids) == {'T2', 'M2', 'ME2'}, cohort.sample_ids
    assert cohort.model2tumor == {
        'M2': 'T2',
        'ME2': 'T2',
    }, cohort.model2tumor
    assert set(cohort.tumor2models.keys()) == {'T2'}, cohort.tumor2models
    assert set(cohort.tumor2models['T2']) == {'M2', 'ME2',}, cohort.tumor2models

def test_cohort_6(links): # exclude M2
    cohort = Cohort(links)
    exclude = ['M2']
    cohort.model_ids = [t for t in cohort.model_ids if t not in exclude]
    assert set(cohort.tumor_ids) == {'T1', 'T2'}, cohort.tumor_ids
    assert set(cohort.model_ids) == {'M1', 'ME2'}, cohort.model_ids
    assert set(cohort.sample_ids) == {'T1', 'T2', 'M1', 'ME2'}, cohort.sample_ids
    assert cohort.model2tumor == {
        'M1': 'T1',
        'ME2': 'T2',
    }, cohort.model2tumor
    assert set(cohort.tumor2models.keys()) == {'T1', 'T2'}, cohort.tumor2models
    assert set(cohort.tumor2models['T1']) == {'M1',}, cohort.tumor2models
    assert set(cohort.tumor2models['T2']) == {'ME2',}, cohort.tumor2models

def test_cohort_7(links): # exclude M2 and ME2
    cohort = Cohort(links)
    exclude = ['M2', 'ME2']
    cohort.model_ids = [t for t in cohort.model_ids if t not in exclude]
    assert set(cohort.tumor_ids) == {'T1'}, cohort.tumor_ids
    assert set(cohort.model_ids) == {'M1'}, cohort.model_ids
    assert set(cohort.sample_ids) == {'T1', 'M1'}, cohort.sample_ids
    assert cohort.model2tumor == {
        'M1': 'T1',
    }, cohort.model2tumor
    assert set(cohort.tumor2models.keys()) == {'T1'}, cohort.tumor2models
    assert set(cohort.tumor2models['T1']) == {'M1',}, cohort.tumor2models

def test_cohort_8(links): # manual edit of mapping blocked
    cohort = Cohort(links)
    try:
        cohort.tumor2models = {} # manual edit
    except AttributeError as e:
        assert 'object has no setter' in str(e)
    else:
        assert False, "AttributeError not raised as expected."

def test_cohort_9(links): # manual edit of mapping blocked
    cohort = Cohort(links)
    try:
        cohort.model2tumor = {} # manual edit
    except AttributeError as e:
        assert 'object has no setter' in str(e)
    else:
        assert False, "AttributeError not raised as expected."

def test_cohort_10(links): # init test
    cohort = Cohort(links)
    assert cohort.sample2sampletype == {
        'T1': 'tumor', 'T2': 'tumor',
        'M1': 'model', 'M2': 'model',
        'ME2': 'model_expanded',
    }, cohort.sample2sampletype

def test_extract_sample_ids_1():
    pair_id = 'HCM-CSHL-0377-C18-85A-01D-A78T-36--HCM-CSHL-0377-C18-10A-01D-A78T-36'
    _tumor_id, _normal_id = 'HCM-CSHL-0377-C18-85A-01D-A78T-36', 'HCM-CSHL-0377-C18-10A-01D-A78T-36'
    _tumor_id_short, _normal_id_short = 'HCM-CSHL-0377-C18-85A', 'HCM-CSHL-0377-C18-10A'
    tumor_id, normal_id, tumor_id_short, normal_id_short = extract_sample_ids(pair_id)
    assert tumor_id == _tumor_id, tumor_id
    assert tumor_id_short == _tumor_id_short, tumor_id_short
    assert normal_id == _normal_id, normal_id
    assert normal_id_short == _normal_id_short, normal_id_short

def test_extract_sample_ids_2():
    pair_id = 'HCM-CSHL-0377-C18-10A-01D-A78T-36--HCM-CSHL-0377-C18-85A-01D-A78T-36'
    _normal_id, _tumor_id = 'HCM-CSHL-0377-C18-10A-01D-A78T-36', 'HCM-CSHL-0377-C18-85A-01D-A78T-36'
    _normal_id_short, _tumor_id_short = 'HCM-CSHL-0377-C18-10A', 'HCM-CSHL-0377-C18-85A'
    tumor_id, normal_id, tumor_id_short, normal_id_short = extract_sample_ids(pair_id)
    assert tumor_id == _tumor_id, tumor_id
    assert tumor_id_short == _tumor_id_short, tumor_id_short
    assert normal_id == _normal_id, normal_id
    assert normal_id_short == _normal_id_short, normal_id_short


def test_interp_nan_1():
    data = np.array([np.nan, np.nan, np.nan, 3.5, np.nan, 
                     np.nan, 6.5,    7.5,    8.5, np.nan])
    expected = np.array([3.5, 3.5, 3.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 8.5])
    result = interp_nan(data)
    assert np.all(expected == result)

def test_interp_nan_2():
    data = np.array([np.nan, np.nan, np.nan])
    expected = np.array([0., 0., 0.,])
    result = interp_nan(data)
    assert np.all(expected == result)

def test_interp_nan_per_chrom_1():
    data = np.array([
        [1, np.nan, 3, 4, np.nan],
        [5, 6, np.nan, 8, 9],
        [4, 3, np.nan, 1, 0],
    ])
    mask = [False, True, True, False, True]
    data = interp_nan_per_chrom(data, mask)
    expected = np.array([[1. , 3. , 3. , 4. , 3. ],
                         [5. , 6. , 7.5, 8. , 9. ],
                         [4. , 3. , 1.5, 1. , 0. ]])
    assert np.all(expected == data)

def test_interp_nan_per_chrom_2():
    data = np.array([
        [1, np.nan, np.nan, np.nan],
        [5, np.nan, np.nan, np.nan],
        [4, np.nan, np.nan, 1],
    ])
    mask = [False, True, True, True]
    data = interp_nan_per_chrom(data, mask)
    expected = np.array([[1., 0., 0., 0.],
                         [5., 0., 0., 0.],
                         [4., 1., 1., 1.]])
    assert np.all(expected == data)

def test_transform_aa_change_1():
    aa_change = 'p.Gly12Val'
    result = transform_aa_change(aa_change)
    assert result == 'p.G12V', result

def test_transform_aa_change_2():
    aa_change = 'p.Gly12Ter'
    result = transform_aa_change(aa_change)
    assert result == 'p.G12*', result

def test_transform_aa_change_3():
    aa_change = 'p.His121ThrfsTer13'
    result = transform_aa_change(aa_change)
    assert result == 'p.H121Tfs*13', result

def test_transform_aa_change_4():
    aa_change = 'p.Phe80='
    result = transform_aa_change(aa_change)
    assert result == 'p.F80=', result