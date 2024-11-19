export init_println
export final_println

@inline function init_println( str ; nhash=9 )
    nstr    = length(str)
    println( "\n" * ("#"^nhash) * " " * str * " " * ("#"^(nhash+11)) )
end

@inline function final_println( str ; nhash=9 )
    nstr    = length(str)
    println( ("#"^nhash) * " " * str * " (finished) " * ("#"^nhash) * "\n" )
end
