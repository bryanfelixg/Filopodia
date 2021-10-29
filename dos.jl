# more initial Conditions
r0[25:75,25:75,1] .= 10;

gr()
data = r0[:,:,1]
heatmap(1:size(data,1),
    1:size(data,2), data,
    c=cgrad([:blue, :white,:red, :yellow]),
    clims =(0,20),
    xlabel="x values", ylabel="y values",
    title="My title")
