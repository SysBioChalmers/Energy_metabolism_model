function [Glc, TCA, OXP, Prd] = calculateProteome(model, prot_cost_info, fluxes, orgName, factor_hy, factor_ly)

if exist('factor_hy', 'var')
    if isempty(factor_hy)
        factor_hy = 1;
    end
else
    factor_hy = 1;
end

if exist('factor_ly', 'var')
    if isempty(factor_ly)
        factor_ly = 1;
    end
else
    factor_ly = 1;
end



if strcmp(orgName, 'Ecoli')
    
    rxn_Glc = {'GLCptspp';'PGI';'PFK';'FBA';'TPI';'GAPD';'PGK';'PGM';'ENO';'PYK'};
    rxn_TCA = {'PDH';'CS';'ACONTa';'ACONTb';'ICDHyr';'AKGDH';'SUCOAS';'FUM';'MDH';'SUCDi';'NADTRHD'};
    rxn_OXP = {'ATPS4rpp';'CYTBO34pp';'NADH16pp'};
    rxn_Prd = {'ACKr';'ACt2rpp';'PTAr'};
    
elseif strcmp(orgName, 'Yeast')

	rxn_Glc = {'GLCt1';'HEX1';'PGI';'PFK';'FBA';'TPI';'GAPD';'PGK';'PGM';'ENO';'PYK';'H2Ot'};
    rxn_TCA = {'PYRt2m';'PDHm';'CSm';'ACONTam';'ACONTbm';'ICDHxm';'AKGDam';'AKGDbm';'GCC2cm';'SUCOASm';'SUCD2u6m';'FUMm';'MDHm'};
    rxn_OXP = {'CYOOm';'CYORu6m';'NADH2u6cm';'NADH2u6m';'ATPS3m';'ATPtmH';'PIt2m'};
    rxn_Prd = {'PYRDC';'ALCD2ir'};
    
end


Glc = 0;
TCA = 0;
OXP = 0;
Prd = 0;

for i = 1:length(model.rxns)
    rxnid = model.rxns{i};
    if contains(rxnid,'_HY')
        id_tmp = rxnid(1:end-3);
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp))/factor_hy;
            flux = fluxes(i);
            
            if ismember(id_tmp,rxn_Glc)
                Glc = Glc + cost * flux;
            elseif ismember(id_tmp,rxn_TCA)
                TCA = TCA + cost * flux;
            elseif ismember(id_tmp,rxn_OXP)
                OXP = OXP + cost * flux;
            elseif ismember(id_tmp,rxn_Prd)
                Prd = Prd + cost * flux;
            end
        end
    elseif contains(rxnid,'_LY')
        id_tmp = rxnid(1:end-3);
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp))/factor_ly;
            flux = fluxes(i);
            
            if ismember(id_tmp,rxn_Glc)
                Glc = Glc + cost * flux;
            elseif ismember(id_tmp,rxn_TCA)
                TCA = TCA + cost * flux;
            elseif ismember(id_tmp,rxn_OXP)
                OXP = OXP + cost * flux;
            elseif ismember(id_tmp,rxn_Prd)
                Prd = Prd + cost * flux;
            end
        end
    elseif contains(rxnid,'_Bio')
        id_tmp = rxnid(1:end-4);
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp));
            flux = fluxes(i);
            
            if ismember(id_tmp,rxn_Glc)
                Glc = Glc + cost * flux;
            elseif ismember(id_tmp,rxn_TCA)
                TCA = TCA + cost * flux;
            elseif ismember(id_tmp,rxn_OXP)
                OXP = OXP + cost * flux;
            elseif ismember(id_tmp,rxn_Prd)
                Prd = Prd + cost * flux;
            end
        end
    end
end

Glc = Glc / 1000;
TCA = TCA / 1000;
OXP = OXP / 1000;
Prd = Prd / 1000;


