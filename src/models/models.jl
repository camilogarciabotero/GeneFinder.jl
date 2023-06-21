const ECOLICDS = begin
    coding = [0.310 0.224 0.199 0.268; 0.251 0.215 0.313 0.221; 0.236 0.308 0.249 0.207; 0.178 0.217 0.338 0.267]
    codinginits = [0.245 0.243 0.273 0.239]
    noncoding = [0.321 0.204 0.200 0.275; 0.282 0.233 0.269 0.215; 0.236 0.305 0.235 0.225; 0.207 0.219 0.259 0.314]
    noncodinginits = [0.262 0.239 0.240 0.259]
    
    TransitionModel(coding, noncoding, codinginits, noncodinginits)
end