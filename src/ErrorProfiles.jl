# TODO READ IN R1 and R2 files, import them into the Dict() consts

struct ErrorRate
    R1::NTuple{150, Float64}
    R2::NTuple{150, Float64}
end

const ErrorRateA = ErrorRate(A_r1_vals, A_r2_vals)


function basequal(rates::ErrorRate, position::UInt16 ,read_number::String)::Tuple{Char, UInt8}
    error = rates[Symbol(read_number)][position]

    

end