# Developers

To develop MicroMag, simply using
```julia
(@v1.7) pkg> dev MicroMag
```
in the Julia REPL. The folder `$JULIA_DEPOT_PATH/dev/MicroMag` should be created and you can modifiy the codes in it. 
For example, we could open the file `src/MicroMag.jl` and add a function `dev_test` 
```julia
function dev_test()
    return "This is a newly added function!"
end
```
We can check it in a new Julia REPL:
```julia
julia> using MicroMag
[ Info: Precompiling MicroMag [8b6b6816-cea2-582c-a99f-83810c20db0f]

julia> MicroMag.dev_test()
"This is a newly added function!"
```

After the modification, we can push our codes into github using 
```bash
git commit -m "we added a dev function" -a
git push
```

