check_base_masks_step2() {
    if [[ ${reads} == 3 ]]; then
        base_mask_SI_SC3Pv2_GEX='Y26n*,I8n*,Y98n*'
        base_mask_SI_SC3Pv2_ADT='Y26n*,I8n*,Y98n*'
        base_mask_SI_SC3Pv3_GEX='Y28n*,I8n*,Y90n*'
        base_mask_SI_SC3Pv3_ADT='Y28n*,I8n*,Y90n*'
        base_mask_SI_SC5P_R2_GEX='Y26n*,I8n*,Y90n*'
        base_mask_SI_SC5P_R2_ADT='Y26n*,I8n*,Y90n*'
        base_mask_SI_SC5P_R2_VDJ='Y26n*,I8n*,Y90n*'
        base_mask_SI_DOGMA_ARCv1_ADT='Y24n*,I8n*,Y90n*'
        base_mask_SI_DOGMA_ARCv1_HTO='Y28n*,I8n*,Y90n*'
    elif [[ ${reads} == 4 ]]; then
        base_mask_SI_SC3Pv2_GEX='Y26n*,I8n*,N*,Y98n*'
        base_mask_SI_SC3Pv2_ADT='Y26n*,I8n*,N*,Y98n*'
        base_mask_SI_SC3Pv3_GEX='Y28n*,I8n*,N*,Y90n*'
        base_mask_SI_SC3Pv3_ADT='Y28n*,I8n*,N*,Y90n*'
        base_mask_SI_SC5P_R2_GEX='Y26n*,I8n*,N*,Y90n*'
        base_mask_SI_SC5P_R2_ADT='Y26n*,I8n*,N*,Y90n*'
        base_mask_SI_SC5P_R2_VDJ='Y26n*,I8n*,N*,Y90n*'
        base_mask_SI_DOGMA_ARCv1_ADT='Y24n*,I8n*,N*,Y90n*'
        base_mask_SI_DOGMA_ARCv1_HTO='Y28n*,I8n*,N*,Y90n*'
        base_mask_DI_SC3Pv2_GEX='Y26n*,I8n*,N*,Y98n*'
        base_mask_DI_SC3Pv2_ADT='Y26n*,I8n*,N*,Y98n*'
        base_mask_DI_SC3Pv3_GEX='Y28n*,I10n*,I10n*,Y90n*'
        base_mask_DI_SC3Pv3_ADT='Y28n*,I10n*,I10n*,Y90n*'
        base_mask_DI_SC3Pv4_GEX='Y28n*,I10n*,I10n*,Y90n*'
        base_mask_DI_SC3Pv4_ADT='Y28n*,I10n*,I10n*,Y90n*'
        base_mask_DI_SC5P_R2_GEX='Y26n*,I10n*,I10n*,Y90n*'
        base_mask_DI_SC5P_R2_ADT='Y26n*,I10n*,I10n*,Y90n*'
        base_mask_DI_SC5P_R2_VDJ='Y26n*,I10n*,I10n*,Y90n*'
        base_mask_DI_SC5P_R2_v3_GEX='Y28n*,I10n*,I10n*,Y90n*'
        base_mask_DI_SC5P_R2_v3_ADT='Y28n*,I10n*,I10n*,Y90n*'
        base_mask_DI_SC5P_R2_v3_VDJ='Y28n*,I10n*,I10n*,Y90n*'
        base_mask_DI_Multiome_ARCv1_GEX='Y28n*,I10n*,I10n*,Y90n*'
        base_mask_DI_Multiome_ARCv1_ATAC='50n*,I8n*,Y24n*,Y49n*'
        base_mask_DI_DOGMA_ARCv1_GEX='Y28n*,I10n*,I10n*,Y90n*'
        base_mask_DI_DOGMA_ARCv1_ATAC='Y100n*,I8n*,Y24n*,Y100n*'
        base_mask_DI_DOGMA_ARCv1_ADT=Y28n*,I8n*,N*,Y90n*
        base_mask_DI_DOGMA_ARCv1_HTO=Y28n*,I8n*,N*,Y90n*
        base_mask_DI_ATAC_ATAC='Y50n*,I8n*,Y16n*,Y50n*'
        base_mask_DI_ASAP_ATAC='Y100n*,I8n*,Y16n*,Y100n*'
        base_mask_DI_ASAP_ADT='Y100n*,I8n*,Y16n*,Y100n*'
        base_mask_DI_ASAP_HTO='Y100n*,I8n*,Y16n*,Y100n*'
        base_mask_DI_ASAP_GENO='Y100n*,I8n*,Y16n*,Y100n*'
    else
        echo -e "\033[0;31mERROR:\033[0m Cannot determine number of reads, check RunInfo.xml file and check_base_masks_step2 criteria"
        exit 1
    fi
}

check_base_masks_step3() {
    local file="$1"
    local run_type="$2"

    # Default values
    local cellranger_command=""
    local index_type=""
    local filter_option=""
    local base_mask=""

    if [[ "${file}" == *_SI_* ]]; then
        index_type='SI'
        filter_option='--filter-single-index'
    elif [[ "${file}" == *_DI_* ]]; then
        index_type='DI'
        filter_option='--filter-dual-index'
    fi

    # Logic for determining the parameters
    if [[ ${file} == CITE_* ]] || [[ ${file} == GEX_* ]]; then
        cellranger_command='cellranger mkfastq'
        if [[ ${file} == *_SI_* ]]; then
            filter_option='--filter-single-index'
            if [[ ${file} == *3prime* ]]; then
                if [[ ${file} == *v2* ]]; then
                    if [[ ${file} == *_GEX ]]; then
                        base_mask=$base_mask_SI_SC3Pv2_GEX
                    elif [[ ${file} == *_ADT ]] || [[ ${file} == *_HTO ]]; then
                        base_mask=$base_mask_SI_SC3Pv2_ADT
                    fi
                elif [[ ${file} == *v3* ]]; then
                    if [[ ${file} == *_GEX ]]; then
                        base_mask=$base_mask_SI_SC3Pv3_GEX
                    elif [[ ${file} == *_ADT ]] || [[ ${file} == *_HTO ]]; then
                        base_mask=$base_mask_SI_SC3Pv3_ADT
                    fi
                fi
            elif [[ ${file} == *SC5P* ]]; then
                if [[ ${file} == *_GEX ]]; then
                    base_mask=$base_mask_SI_SC5P_R2_GEX
                elif [[ ${file} == *_ADT ]] || [[ ${file} == *_HTO ]]; then
                    base_mask=$base_mask_SI_SC5P_R2_ADT
                elif [[ ${file} == *_VDJ* ]]; then
                    base_mask=$base_mask_SI_SC5P_R2_VDJ
                fi
            fi
        elif [[ ${file} == *_DI_* ]]; then
            filter_option='--filter-dual-index'
            if [[ ${file} == *3prime* ]]; then
                if [[ ${file} == *v2* ]]; then
                    if [[ ${file} == *_GEX ]]; then
                        base_mask=$base_mask_DI_SC3Pv2_GEX
                    elif [[ ${file} == *_ADT ]] || [[ ${file} == *_HTO ]]; then
                        base_mask=$base_mask_DI_SC3Pv2_ADT
                    fi
                elif [[ ${file} == *v3* ]]; then
                    if [[ ${file} == *_GEX ]]; then
                        base_mask=$base_mask_DI_SC3Pv3_GEX
                    elif [[ ${file} == *_ADT ]] || [[ ${file} == *_HTO ]]; then
                        base_mask=$base_mask_DI_SC3Pv3_ADT
                    fi
                elif [[ ${file} == *v4* ]]; then
                    if [[ ${file} == *_GEX ]]; then
                        base_mask=$base_mask_DI_SC3Pv4_GEX
                    elif [[ ${file} == *_ADT ]] || [[ ${file} == *_HTO ]]; then
                        base_mask=$base_mask_DI_SC3Pv4_ADT
                    fi
                fi
            elif [[ ${file} == *SC5P* ]]; then
                if [[ ${file} == *v2* ]]; then
                    if [[ ${file} == *_GEX ]]; then
                        base_mask=$base_mask_DI_SC5P_R2_GEX
                    elif [[ ${file} == *ADT ]] || [[ ${file} == *_HTO ]]; then
                        base_mask=$base_mask_DI_SC5P_R2_ADT
                    elif [[ ${file} == *_VDJ* ]]; then
                        base_mask=$base_mask_DI_SC5P_R2_VDJ
                    fi
                elif [[ ${file} == *v3* ]]; then
                    if [[ ${file} == *_GEX ]]; then
                        base_mask=$base_mask_DI_SC5P_R2_v3_GEX
                    elif [[ ${file} == *_ADT ]] || [[ ${file} == *_HTO ]]; then
                        base_mask=$base_mask_DI_SC5P_R2_v3_ADT
                    elif [[ ${file} == *_VDJ* ]]; then
                        base_mask=$base_mask_DI_SC5P_R2_v3_VDJ
                    fi
                fi
            fi
        fi
    elif [[ ${file} == Multiome_* ]]; then
        cellranger_command='cellranger mkfastq'
        filter_option='--filter-dual-index'
        if [[ ${file} == *_GEX ]]; then
            base_mask=$base_mask_DI_Multiome_ARCv1_GEX
        elif [[ ${file} == *_ATAC ]]; then
            cellranger_command='cellranger-atac mkfastq'
            base_mask=$base_mask_DI_Multiome_ARCv1_ATAC
        elif [[ ${file} == *_ADT ]]; then
            base_mask=$base_mask_DI_Multiome_ARCv1_ADT
        elif [[ ${file} == *_HTO ]]; then
            base_mask=$base_mask_DI_Multiome_ARCv1_HTO
        fi
    elif [[ ${file} == DOGMA_* ]]; then
        if [[ ${file} == *_SI_* ]]; then
            cellranger_command='cellranger mkfastq'
            filter_option='--filter-single-index'
            if [[ ${file} == *_ADT ]]; then
                base_mask=$base_mask_SI_DOGMA_ARCv1_ADT
            elif [[ ${file} == *_HTO ]]; then
                base_mask=$base_mask_SI_DOGMA_ARCv1_HTO
            fi
        elif [[ ${file} == *_DI_* ]]; then
            cellranger_command='cellranger mkfastq'
            filter_option='--filter-dual-index'
            if [[ ${file} == *_GEX ]]; then
                base_mask=$base_mask_DI_DOGMA_ARCv1_GEX
            elif [[ ${file} == *_ATAC ]]; then
                cellranger_command='cellranger-atac mkfastq'
                base_mask=$base_mask_DI_DOGMA_ARCv1_ATAC
            elif [[ ${file} == *_ADT ]]; then
                base_mask=$base_mask_DI_DOGMA_ARCv1_ADT
            elif [[ ${file} == *_HTO ]]; then
                base_mask=$base_mask_DI_DOGMA_ARCv1_HTO
            fi
        fi
    elif [[ ${file} == ASAP_* ]]; then
        cellranger_command='cellranger-atac mkfastq'
        filter_option='--filter-dual-index'
        if [[ ${file} == *_ATAC ]]; then
            base_mask=$base_mask_DI_ASAP_ATAC
        elif [[ ${file} == *_ADT ]]; then
            base_mask=$base_mask_DI_ASAP_ADT
        elif [[ ${file} == *_HTO ]]; then
            base_mask=$base_mask_DI_ASAP_HTO
        elif [[ ${file} == *_GENO ]]; then
            base_mask=$base_mask_DI_ASAP_GENO
        fi
    elif [[ ${file} == ATAC_* ]]; then
        cellranger_command='cellranger-atac mkfastq'
        filter_option='--filter-dual-index'
        base_mask=$base_mask_DI_ATAC_ATAC
    else
        echo -e "\033[0;31mERROR:\033[0m Cannot determine base mask for ${file}, please check path"
        exit 1
    fi

    # Export the variables if needed
    echo "${cellranger_command// /.}" "${index_type// /.}" "${filter_option// /.}" "${base_mask// /.}"
}