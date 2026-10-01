include { chemistry_field } from './chemistry'

// CITE/GEX masks depend only on index type and chemistry family; VDJ needs a 5' chemistry.
def get_override_cycles(Map sp, int num_reads) {
    def masks = [
        'SI_SC3Pv2':              [3: 'Y26N*;I8N*;Y98N*',  4: 'Y26N*;I8N*;N*;Y98N*'],
        'SI_SC3Pv3':              [3: 'Y28N*;I8N*;Y90N*',  4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC3Pv4':              [3: 'Y28N*;I8N*;Y90N*',  4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC5P':                [3: 'Y26N*;I8N*;Y90N*',  4: 'Y26N*;I8N*;N*;Y90N*'],
        'DI_SC3Pv2':              [4: 'Y26N*;I8N*;N*;Y98N*'],
        'DI_SC3Pv3':              [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC3Pv4':              [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC5P':                [4: 'Y26N*;I10N*;I10N*;Y90N*'],
        'DI_SC5Pv3':              [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'SI_DOGMA_ARCv1_ADT':     [3: 'Y24N*;I8N*;Y90N*',  4: 'Y24N*;I8N*;N*;Y90N*'],
        'SI_DOGMA_ARCv1_HTO':     [3: 'Y28N*;I8N*;Y90N*',  4: 'Y28N*;I8N*;N*;Y90N*'],
        'DI_Multiome_ARCv1_GEX':  [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_Multiome_ARCv1_ATAC': [4: 'Y50N*;I8N*;Y24N*;Y49N*'],
        'DI_DOGMA_ARCv1_GEX':     [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_DOGMA_ARCv1_ATAC':    [4: 'Y100N*;I8N*;Y24N*;Y100N*'],
        'DI_DOGMA_ARCv1_ADT':     [4: 'Y28N*;I8N*;N*;Y90N*'],
        'DI_DOGMA_ARCv1_HTO':     [4: 'Y28N*;I8N*;N*;Y90N*'],
        'DI_ATAC_ATAC':           [4: 'Y50N*;I8N*;Y16N*;Y50N*'],
        'DI_ASAP_ATAC':           [4: 'Y100N*;I8N*;Y16N*;Y100N*'],
        'DI_ASAP_ADT':            [4: 'Y100N*;I8N*;Y16N*;Y100N*'],
        'DI_ASAP_HTO':            [4: 'Y100N*;I8N*;Y16N*;Y100N*'],
        'DI_Flex-v2_GEX':         [4: 'Y*;I10;I10;Y*'],
    ]

    def mod = sp.modality.replaceAll(/^VDJ-[TB]$/, 'VDJ').replaceAll(/^CRISPR$/, 'GEX')
    def family = sp.assay in ['CITE', 'GEX'] ? chemistry_field(sp.chemistry, 'family') : null
    def key = family && (mod in ['GEX', 'ADT', 'HTO'] || (mod == 'VDJ' && family.startsWith('SC5P'))) ? "${sp.index_type}_${family}"
            : sp.assay in ['CITE', 'GEX'] ? null
            : sp.assay == 'Flex'          ? 'DI_Flex-v2_GEX'
            : sp.assay == 'Multiome'      ? "DI_Multiome_ARCv1_${mod}"
            : sp.assay == 'DOGMA'         ? (sp.modality == 'ATAC' ? 'DI_DOGMA_ARCv1_ATAC' : "${sp.index_type}_DOGMA_ARCv1_${mod}")
            : sp.assay == 'ASAP'          ? "DI_ASAP_${mod}"
            : sp.assay == 'ATAC'          ? 'DI_ATAC_ATAC'
            : null
    def oc = masks[key?.toString()]
    if (!oc)
        error "Cannot determine OverrideCycles: assay=${sp.assay} chem=${sp.chemistry} index_type=${sp.index_type} modality=${sp.modality} (key=${key})"
    oc = oc[num_reads]
    if (!oc) return null

    if (sp.index_len == 8) oc = oc.replaceAll(/I10N\*/, 'I8N2')
    if (num_reads == 4 && !sp.is_dual) oc = apply_si_on_di_correction(oc, sp.index_len)
    return oc
}

def apply_si_on_di_correction(String oc, Integer seq_len, Integer index1_cycles = 10) {
    def parts     = oc.split(';') as List
    def i_indices = (0..<parts.size()).findAll { i -> parts[i].startsWith('I') }

    if (i_indices.isEmpty()) return oc

    def fixed = (0..<parts.size()).collect { idx ->
        if (idx == i_indices[0]) {
            def remaining = index1_cycles - seq_len
            remaining > 0 ? "I${seq_len}N${remaining}" : "I${seq_len}"
        } else if (i_indices.size() > 1 && idx == i_indices[1]) {
            'N*'
        } else {
            parts[idx]
        }
    }

    return fixed.join(';')
}
