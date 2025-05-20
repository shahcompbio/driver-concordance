import numpy as np
import pandas as pd
import anndata as ad

from gdan_hcmi_tools.copynumber import get_cnv_events_hcmi, compute_cn_distance, compute_mean_distance, eval_cn_event, gene_over_min_cn, has_hlamp, rebin

def test_rebin():
    data = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1', 'chr2', 'chr2', 'chr3'],
        'Start.bp': [0, 100, 0, 100, 0],
        'End.bp': [100, 200, 100, 200, 100],
        'value': [1, 2, 3, 4, 5]
    })
    data['Chromosome'] = data['Chromosome'].str.replace('chr', '')
    data = data[['Chromosome', 'Start.bp', 'End.bp', 'value']]
    bin_size = 50
    cols = ['value']
    rebinned = rebin(data, bin_size, cols)
    # expected = pd.DataFrame({
    #     'chromosome': ['1', '1', '1', '1', '2', '2', '2', '2', '3', '3'],
    #     'start': [0, 50, 100, 150] + [0, 50, 100, 150] + [0, 50],
    #     'end': [50, 100, 150, 200] + [50, 100, 150, 200] + [50, 100],
    #     'value': [1.02, 1.02, 2.04, 2.04] + [3.06, 3.06, 4.08, 4.08] + [5.10, 5.10],
    # })
    # write `expected`, based on the AssertionError output above
    expected = pd.DataFrame({
        'chromosome': ['1', '1', '1', '1', '2', '2', '2', '2', '3', '3'],
        'start': [0, 50, 100, 150] + [0, 50, 100, 150] + [0, 50],
        'end': [50, 100, 150, 200] + [50, 100, 150, 200] + [50, 100],
        'value': [1.02, 1.02, 2.04, 2.04] + [3.06, 3.06, 4.08, 4.08] + [5.10, 5.10],
    })
    assert np.allclose(expected.values[:, 1:].astype(float), rebinned.values[:, 1:].astype(float), atol=0.1), rebinned


def make_gene_level_anndata_for_pytest(samples, genes, layer_data,
                                       obs_data=[['tumor', 'PAAD', 1]], layer_name='total_cn'):
    if type(layer_data) == list: 
        layer_data = np.array(layer_data)
    assert len(obs_data) == len(samples), (samples, obs_data)
    obs = pd.DataFrame(index=samples, columns=['sample_type', 'cancer_type', 'ploidy'], 
                       data=obs_data)
    var = pd.DataFrame(index=genes)
    gdata = ad.AnnData(layers={layer_name:layer_data}, obs=obs, var=var)
    return gdata

def test_gene_over_min_cn_1():
    assert gene_over_min_cn(gene_cn=4, gene_cn_remixt=4, min_gene_cn=2)

def test_gene_over_min_cn_2():
    assert gene_over_min_cn(gene_cn=1, gene_cn_remixt=4, min_gene_cn=2)

def test_has_hlamp_1():
    assert has_hlamp(gene_cn=6, hlamp_cutoff=5, ploidy=2, call_cn_amp_cn=6, min_ploidy_cutoff=2, debug=False)

def test_eval_cn_event_1():
    event = eval_cn_event(gene_cn=6., gene_cn_remixt=6.,
                          ploidy_consensus=2., ploidy_remixt=2.,
                          hlamp_cutoff_consensus=4., hlamp_cutoff_remixt=4.,
                          min_gene_cn=2., call_cn_amp_cn=6., min_ploidy_cutoff=2.)
    assert event == 'hlamp', event

def test_eval_cn_event_2():
    event = eval_cn_event(gene_cn=2., gene_cn_remixt=6.,
                          ploidy_consensus=2., ploidy_remixt=2.,
                          hlamp_cutoff_consensus=4., hlamp_cutoff_remixt=4.,
                          min_gene_cn=2., call_cn_amp_cn=6., min_ploidy_cutoff=2.)
    assert event == 'hlamp', event

def test_eval_cn_event_3():
    event = eval_cn_event(gene_cn=0.3, gene_cn_remixt=0.3,
                          ploidy_consensus=2., ploidy_remixt=2.,
                          hlamp_cutoff_consensus=4., hlamp_cutoff_remixt=4.,
                          homdel_cutoff=0.5)
    assert event == 'homdel', event

def test_eval_cn_event_4():
    event = eval_cn_event(gene_cn=3, gene_cn_remixt=0.3,
                          ploidy_consensus=2., ploidy_remixt=2.,
                          hlamp_cutoff_consensus=4., hlamp_cutoff_remixt=4.,
                          homdel_cutoff=0.5)
    assert event == 'homdel', event

def test_eval_cn_event_5():
    event = eval_cn_event(gene_cn=6, gene_cn_remixt=6,
                          ploidy_consensus=3., ploidy_remixt=3.,
                          hlamp_cutoff_consensus=6., hlamp_cutoff_remixt=6.,
                          homdel_cutoff=0.5)
    assert event == 'hlamp', event

# test for a case the gene_cn is below the min_gene_cn
def test_eval_cn_event_6():
    event = eval_cn_event(gene_cn=1.5, gene_cn_remixt=1.5,
                          ploidy_consensus=2., ploidy_remixt=2.,
                          hlamp_cutoff_consensus=4., hlamp_cutoff_remixt=4.,
                          min_gene_cn=2., call_cn_amp_cn=6., min_ploidy_cutoff=2.)
    assert event is None, event

def test_get_cnv_events_hcmi_1():
    gene_sets = {'hlamp':{'A1CF'}, 'homdel':{}}
    obs_data = [['tumor', 'PAAD', 2.1]]
    gdata = make_gene_level_anndata_for_pytest(['sample'], ['A1CF'], [[2*2.1+0.5+0.1]], obs_data=obs_data)
    rgdata = gdata.copy()
    cnv_events = get_cnv_events_hcmi(gdata, rgdata, gene_sets, b_consensus=0.)
    expected_amp = [[1]]
    assert np.all(cnv_events['hlamp'] == expected_amp)

def test_get_cnv_events_hcmi_2():
    gene_sets = {'hlamp':{'A1CF'}, 'homdel':{}}
    obs_data = [['tumor', 'PAAD', 2.5]]
    gdata = make_gene_level_anndata_for_pytest(['sample'], ['A1CF'], [[4.]], obs_data=obs_data)
    rgdata = gdata.copy()
    cnv_events = get_cnv_events_hcmi(gdata, rgdata, gene_sets)
    expected_amp = [[0]]
    assert np.all(cnv_events['hlamp'] == expected_amp), cnv_events['hlamp']

def test_get_cnv_events_hcmi_3():
    gene_sets = {'homdel':{'CDKN2A'}, 'hlamp':{}}
    obs_data = [['tumor', 'GBM', 2.1]]
    gdata = make_gene_level_anndata_for_pytest(['sample'], ['CDKN2A'], [[0]], obs_data=obs_data)
    rgdata = make_gene_level_anndata_for_pytest(['sample'], ['CDKN2A'], [[1]], obs_data=obs_data)
    cnv_events = get_cnv_events_hcmi(gdata, rgdata, gene_sets)
    expected_del = [[1]]
    assert np.all(cnv_events['homdel'] == expected_del), cnv_events['homdel']

def test_get_cnv_events_hcmi_4():
    gene_sets = {'homdel':{'CDKN2A'}, 'hlamp':{}}
    obs_data = [['tumor', 'GBM', 2.1]]
    gdata = make_gene_level_anndata_for_pytest(['sample'], ['CDKN2A'], [[1]], obs_data=obs_data)
    rgdata = make_gene_level_anndata_for_pytest(['sample'], ['CDKN2A'], [[0]], obs_data=obs_data)
    cnv_events = get_cnv_events_hcmi(gdata, rgdata, gene_sets)
    expected_del = [[1]]
    assert np.all(cnv_events['homdel'] == expected_del), cnv_events['homdel']

def test_compute_cn_distance():
    ploidies = [0, 0]
    paired_data = pd.DataFrame({
        'start': [0, 100, 200],
        'end': [100, 200, 1000],
        'major_raw_e_tumor': [3, 3, 2.4],
        'minor_raw_e_tumor': [2, 2, 2],
        'major_raw_e_model': [1, 1, 1],
        'minor_raw_e_model': [1, 1, 1],
    })
    paired_data, logs = compute_cn_distance(paired_data, ploidies)
    paired_data['major_distance']
    major_distance = compute_mean_distance(paired_data, metric='major_distance')
    minor_distance = compute_mean_distance(paired_data, metric='minor_distance')
    assert paired_data['major_distance'].round(1).tolist() == [1.0, 1.0, 0.4]
    assert paired_data['minor_distance'].tolist() == [0, 0, 0]
    assert major_distance == (800 * 0.4 + 100*1 + 100*1) / (800+100+100), major_distance
    assert minor_distance == 0, minor_distance