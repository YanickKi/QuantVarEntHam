using Yao 

function whole()
    A = [put(2,1=>X), put(2,2=>X)]
    B = mat.(A)
    println(typeof(B))

end 

whole()