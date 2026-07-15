using CSV
using DataFrames

"""
    convert_txt_to_csv(filename::String)

Checks if a CSV version of `filename` exists in the data directory. 
If not, reads the TXT version and converts it to CSV without overwriting.
This is useful for Honda tables which are split using spaces. 
"""
function convert_txt_to_csv(filename::String;data_dir= joinpath(@__DIR__, "..", "data"))
    
    # construct paths

    txt_path = joinpath(data_dir, "$(filename).txt")
    csv_path = joinpath(data_dir, "$(filename).csv")

    # check if the CSV file already exists to avoid overwriting it
    if isfile(csv_path)
        println("Skipped: '$(filename).csv' already exists. Keeping the original version.")
    end

    # check if the source TXT file actually exists before reading
    if !isfile(txt_path)
        error("Error: Source file '$(filename).txt' not found in $(data_dir)")
    end

    # read the TXT file and write the new CSV
    df = CSV.read(txt_path, DataFrame, delim=' ', ignorerepeated=true)
    CSV.write(csv_path, df)
    println("Successfully created: '$(filename).csv'")

    return df

end

#convert_txt_to_csv("nuflux")