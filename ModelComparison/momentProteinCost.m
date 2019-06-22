function prot_cost_info = momentProteinCost(model)

k_eff = 65; % /s

prot_cost_info = struct();
prot_cost_info.id = model.rxns;
prot_cost_info.value = zeros(length(model.rxns),1);

[num, txt, ~] = xlsread('MW_ecoli.xlsx');
MW_gene = txt(2:end,1);
MW_value = num;

for i = 1:length(model.rxns)
    z = model.rules{i}; % get gpr
    if isempty(z)%if is spontaneous reaction
        prot_cost_info.value(i) = 0;
    else%if is enzymatic reaction
        if ~contains(z,' | ') && ~contains(z,' & ')%if has a single enzyme
            mw = MW_value(ismember(MW_gene,z));
            prot_cost_info.value(i) = mw/k_eff/3600;
        else%if has multiple genes, either isozymes or complexes
            z = z(2:length(z)-1);%delete parentheses
            z = strtrim(z);%delete leading and trailing whitespace if it has
            if ~contains(z,' | ')%if only has complex
                z = strsplit(z,' & ');%split isozymes
                mw = sum(MW_value(ismember(MW_gene,z)));
                prot_cost_info.value(i) = mw/k_eff/3600;
            else%if has isozymes, and has or does not have complex
                if ~contains(z,' & ')%if only has isozymes
                    z = strsplit(z,' | ');%split isozymes
                    mw = min(MW_value(ismember(MW_gene,z)));
                    prot_cost_info.value(i) = mw/k_eff/3600;
                else%if has both isozymes and complex
                    if contains(z,' | ')
                    % if GPR is made up of iso-complexes, e.g., '( A & B ) | ( C & D )'
                        z = strsplit(z,' | ');%split isozymes
                        mw_list = zeros(1,length(z));
                        for j = 1:length(z)
                            z_tmp = z{j};
                            if contains(z_tmp,' & ')
                            	z_tmp = z_tmp(2:length(z_tmp)-1);%delete parentheses
                                z_tmp = strtrim(z_tmp);%delete leading and trailing whitespace if it has
                                z_tmp_tmp = strsplit(z_tmp,' & ');%split isozymes
                                mw_list(1,j) = sum(MW_value(ismember(MW_gene,z_tmp_tmp)));
                            else
                                mw_list(1,j) = MW_value(ismember(MW_gene,z_tmp));
                            end
                            mw = min(mw_list);
                            prot_cost_info.value(i) = mw/k_eff/3600;
                        end
                    end
                end
            end
        end
    end
end
