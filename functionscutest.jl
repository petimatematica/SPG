using CUTEst


#custom_filter = x -> x["origin"] == "real"
#custom_filter = x -> x["origin"] == "academic"
#custom_filter = x -> x["origin"] == "modelling"
problems1 = CUTEst.select(objtype="quadratic", contype="bounds")#, custom_filter=custom_filter)
problems2 = CUTEst.select(objtype = "sum_of_squares", contype = "bounds")#, custom_filter=custom_filter)

problems = [problems1, problems2]
length(problems)
println(problems)