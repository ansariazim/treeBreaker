%% modify a phylogeny randomly.

% first of all let's read the tree.
sim_tree = phytreeread('tree1000.nwk');
% get the data in a useful format.
sim_tree_data = get(sim_tree);
% get pairwise distances
sim_tree_pdist = pdist(sim_tree);
% get the data in format that can be plotted by dendrogram function.
sim_tree_Z = linkage(sim_tree_pdist');
figure
dendrogram(sim_tree_Z,0);

% % % my_distsMat = squareform(my_pdist);


%%%%% get the new tree.

% change the distances such that new_dist = unif(0.5*old_distance,1.5*old_distance)
new_pdist = sim_tree_pdist;
for i=1:numel(new_pdist)
    new_pdist(i) = unifrnd(0.5*sim_tree_pdist(i),1.5*sim_tree_pdist(i));
end

% Get the UPGMA tree from the pairwise distances.
new_tree_Z = linkage(new_pdist');
figure
dendrogram(new_tree_Z,0)
% convert it into phytree format.
new_tree = phytree(new_tree_Z,get(sim_tree,'LEAFNAMES'));
% get the newick str.
new_tree_newick = getnewickstr(new_tree);
% get the tree data in a useful format.
new_tree_data = get(new_tree);

% write the newick formatted tree.
fid = fopen('tree1000_mod.nwk','w');
fprintf(fid,'%s',new_tree_newick);
fclose(fid);
