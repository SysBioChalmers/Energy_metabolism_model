%% reformulateEcoliModel 
function new_model = reformulateEcoliModel(model)

% Generate a matrix for splitting reversible reactions.
matrix2 = struct();
matrix2.RxnList = cell(0,1);
matrix2.CompoList = cell(0,1);
matrix2.CoeffList = zeros(0,1);
matrix2.CatalystList = cell(0,1);
matrix2.LBList = zeros(0,1);
matrix2.UBList = zeros(0,1);
matrix2.protCost = zeros(0,1);

for i = 1:length(model.rxns)

    idx_substrate = model.S(:,i) < 0;
    idx_product = model.S(:,i) > 0;
    CompoS = model.mets(idx_substrate);
    CompoP = model.mets(idx_product);
    CoeffS = model.S(idx_substrate,i);
    CoeffP = model.S(idx_product,i);
    CompoList = [CompoS;CompoP];
    CoeffList = [CoeffS;CoeffP];
    n = length(CompoList);
    Rxn = model.rxns(i);
    Catalyst = model.rules(i);
    lb_tmp = model.lb(i);
    ub_tmp = model.ub(i);
    protcost_tmp = model.protCost(i);
    RxnList = repmat(Rxn,n,1);
    CatalystList = repmat(Catalyst,n,1);
    LBList = repmat(lb_tmp,n,1);
    UBList = repmat(ub_tmp,n,1);
    ProtCostList = repmat(protcost_tmp,n,1);
    
    x = length(matrix2.RxnList);%count rows in matrix2
    if isempty(CatalystList{1})%if is spontaneous reaction
        matrix2.RxnList(x+1:x+n,1) = RxnList;
        matrix2.CompoList(x+1:x+n,1) = CompoList;
        matrix2.CoeffList(x+1:x+n,1) = CoeffList;
        matrix2.CatalystList(x+1:x+n,1) = CatalystList;
        matrix2.LBList(x+1:x+n,1) = LBList;
        matrix2.UBList(x+1:x+n,1) = UBList;
        matrix2.protCost(x+1:x+n,1) = ProtCostList;
    else%if is enzymatic reaction
        if LBList(1) >= 0%if is not reversible
            matrix2.RxnList(x+1:x+n,1) = RxnList;
            matrix2.CompoList(x+1:x+n,1) = CompoList;
            matrix2.CoeffList(x+1:x+n,1) = CoeffList;
            matrix2.CatalystList(x+1:x+n,1) = CatalystList;
            matrix2.LBList(x+1:x+n,1) = LBList;
            matrix2.UBList(x+1:x+n,1) = UBList;
            matrix2.protCost(x+1:x+n,1) = ProtCostList;
        else%if is reversible
            %add forward rows
            if UBList(1) > 0
                rxnname_tmp = strcat(RxnList{1},'_fwd');
                matrix2.RxnList(x+1:x+n,1) = repmat({rxnname_tmp},n,1);
                matrix2.CompoList(x+1:x+n,1) = CompoList;
                matrix2.CoeffList(x+1:x+n,1) = CoeffList;
                matrix2.CatalystList(x+1:x+n,1) = CatalystList;
                matrix2.LBList(x+1:x+n,1) = zeros(n,1);
                matrix2.UBList(x+1:x+n,1) = UBList;
                matrix2.protCost(x+1:x+n,1) = ProtCostList;
            end
            x = length(matrix2.RxnList);
            %add reverse rows
            rxnname_tmp = strcat(RxnList{1},'_rvs');
            matrix2.RxnList(x+1:x+n,1) = repmat({rxnname_tmp},n,1);
            matrix2.CompoList(x+1:x+n,1) = CompoList;
            matrix2.CoeffList(x+1:x+n,1) = -1*CoeffList;
            matrix2.CatalystList(x+1:x+n,1) = CatalystList;
            matrix2.LBList(x+1:x+n,1) = zeros(n,1);
            matrix2.UBList(x+1:x+n,1) = -1*LBList;
            matrix2.protCost(x+1:x+n,1) = ProtCostList;
        end
    end
end

matrixSplit = matrix2;

% Converted to COBRA model
new_model = struct();
new_model.rxns = cell(0,1);
new_model.lb = zeros(0,1);
new_model.ub = zeros(0,1);
new_model.mets = cell(0,1);
new_model.S = sparse(0,0);
new_model.b = zeros(0,1);
new_model.rxnGeneMat = sparse(0,0);
new_model.c = zeros(0,1);
new_model.rules = cell(0,1);
new_model.genes = cell(0,1);
new_model.csense = char();

new_model.metFormulas = cell(0,1);
new_model.metNames = cell(0,1);
new_model.metKEGGID = cell(0,1);
new_model.metChEBIID = cell(0,1);
new_model.rxnNames = cell(0,1);

new_model.osenseStr = model.osenseStr;

if isfield(model,'compNames')
    new_model.compNames = model.compNames;
end
if isfield(model,'comps')
    new_model.comps = model.comps;
end

UnqRxnList = unique(matrixSplit.RxnList);
for i = 1:length(UnqRxnList)
    
    id_tmp = UnqRxnList(i);
    idx = ismember(matrixSplit.RxnList,id_tmp);
    
    RxnList = matrixSplit.RxnList(idx);
    CompoList = matrixSplit.CompoList(idx);
    CoeffList = matrixSplit.CoeffList(idx);
    CatalystList = matrixSplit.CatalystList(idx);
    LBList = matrixSplit.LBList(idx);
    UBList = matrixSplit.UBList(idx);
    protCostList = matrixSplit.protCost(idx);
    
    rxn_tmp = RxnList{1};
    catalyst_tmp = CatalystList{1};
    lb_tmp = unique(LBList);
    ub_tmp = unique(UBList);
    protcost_tmp = unique(protCostList);
    protcost_tmp_str = num2str(protcost_tmp);
    
    % add metabolites
    [isInModelidx,~] = ismember(CompoList,new_model.mets);
    NotInModel = CompoList(~isInModelidx);
    idx_met_tmp = ismember(model.mets,NotInModel);
	metid_tmp = model.mets(idx_met_tmp);
	metname_tmp = model.metNames(idx_met_tmp);
	metformula_tmp = model.metFormulas(idx_met_tmp);
	ChEBIiD_tmp = model.metChEBIID(idx_met_tmp);
	KEGGid_tmp = model.metKEGGID(idx_met_tmp);
    new_model = addMetabolite(new_model,metid_tmp,...
                              'metName',metname_tmp,...
                              'metFormula',metformula_tmp,...
                              'ChEBIID',ChEBIiD_tmp,...
                              'KEGGId',KEGGid_tmp);
    
    % add reactions
    new_model = addReaction(new_model,rxn_tmp,...
                            'metaboliteList',CompoList,...
                            'stoichCoeffList',CoeffList,...
                            'lowerBound',lb_tmp,...
                            'upperBound',ub_tmp,...
                            'geneRule',catalyst_tmp,...
                            'subSystem',protcost_tmp_str);
end

new_model.protCost = cellfun(@(x) str2double(x),...
                            new_model.subSystems,'UniformOutput',false);
new_model.protCost = cell2mat(new_model.protCost);                        