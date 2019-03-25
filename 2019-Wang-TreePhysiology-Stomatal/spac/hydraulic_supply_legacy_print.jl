# print the lagcies of root, stem, and leaf
function TreeLegaciesPrint(legacies)
    print_with_color(:red, "\nThe legacies in root, stem, and leaf are listed below:\n")
    println("RootTension\tRootKratio\tStemTension\tStemKratio\tLeafTension\tLeafKratio")
    for i in 1:200
        println(legacies[1,i,1],"\t",legacies[1,i,2],"\t",legacies[2,i,1],"\t",legacies[2,i,2],"\t",legacies[3,i,1],"\t",legacies[3,i,2])
    end
    print_with_color(:red, "Finished printing the legacies!\n")
end
