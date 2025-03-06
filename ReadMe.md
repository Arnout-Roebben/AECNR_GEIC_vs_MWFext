# A comparative analysis of generalised echo and interference cancelling and extended multichannel Wiener filtering for combined noise reduction and acoustic echo cancellation
# License
This work is licensed under the [MIT LICENSE](LICENSE.md). By downloading and/or installing this software and associated files on your computing system you agree to use the software under the terms and conditions as specified in the license agreement.

If this code has been useful for you, please cite [[1]](#References).

# About
This repository [[2]](#References) contains the MATLAB code associated with [[1]](#References). For the problem of combined noise reduction (NR) and acoustic echo cancellation (AEC), where a comparative analysis is conducted between the generalised echo and interference canceller (GEIC) [[4]](#References) and the extended multichannel Wiener filter (MWF<sub>ext</sub>) [[5]](#References). To this end, in [[1]](#References) 1) general nonlinear echo paths, and 2) voice activity detector (VAD) error effects are considered in a comparative analysis between the GEIC and MWF<sub>ext</sub>.

This repository provides example code for the GEIC and MWF<sub>ext</sub>. 

The code has been developed and tested in MATLAB R2024a.

The manuscript of [[1]](#References) has also been released as a preprint [[3]](#References).

# File structure
* [Main.m](Main.m): Main file to run the code.
* [LICENSE](LICENSE.md): License file.
* [ReadMe.md](ReadMe.md): ReadMe file.
* [Audio](Audio): Folder containing the audio files.
    - [impulse.mat](Audio/impulse.mat): File containing the impulse responses between the desired speech source and the microphones.
    - [sig.mat](Audio/sig.mat): File containing the desired speech, near-end room noise, echo, and loudspeaker signals. The speech and speech component in the echo are taken from the [VCTK corpus](https://datashare.ed.ac.uk/handle/10283/3443), which is made available under the Open Data Commons Attribution License (ODC-By) v1.0.
* [Util](Util): Auxiliary code.
    - [AEC+NR](Util/AEC+NR): Folder containing the auxiliary code for combined AEC and NR.  
        + [process.m](Util/AEC+NR/process.m): Processing using the GEIC or MWF<sub>ext</sub>.
        + [process_GEIC.m](Util/AEC+NR/process_GEIC.m): Processing using the GEIC with ground-truth relative transfer functions available.        
        + [process_GEIC_GEVD.m](Util/AEC+NR/process_GEIC_GEVD.m): Processing using the GEIC with estimated relative transfer functions using the covariance whitening method [[6]](#References).       
        + [process_MWFext.m](Util/AEC+NR/process_MWFext.m): Processing using the MWF<sub>ext</sub>.   
    - [Filter](Util/Filter): Folder containing the auxiliary code for filter application.
        + [applyFilterMultichannel.m](Util/Filter/applyFilterMultichannel.m): Applies a multichannel filer.
    - [Frequency_transformation](Util/Frequency_transformation): Folder containing the auxiliary code for the conversion between time- and frequency-domain.
        + [WOLA_analysis.m](Util/Frequency_transformation/WOLA_analysis.m): Weighted overlapp add (WOLA) analysis filterbank.
        + [WOLA_synthesis.m](Util/Frequency_transformation/WOLA_synthesis.m): WOLA synthesis filterbank.  
    - [GEVD](Util/GEVD): Folder containing the auxiliary code for the generalised eigenvalue decomposition (GEVD).
        + [updateDifferenceCorrelation.m](Util/GEVD/updateDifferenceCorrelation.m): Computes the difference between two correlation matrices using a GEVD approximation.
    - [Metrics](Util/Metrics): Folder containing the auxiliary code for the evaluation of the algorithms.
        + [compute_metrics.m](Util/Metrics/compute_metrics.m): Computes the signal-to-noise ratio (SNR), echo return 
        loss enhancement (ERLE), and speech distortion (SD).
        + [ERLE.m](Util/Metrics/ERLE.m): Computes the ERLE.
        + [SD.m](Util/Metrics/SD.m): Computes the SD.
        + [SNR.m](Util/Metrics/SNR.m): Computes the SNR.
    - [VAD](Util/VAD): Folder containing the auxiliary code for the voice activity detection (VAD).
        + [VAD.m](Util/VAD/VAD.m): Computes the VAD in the STFT domain.        

# References
[1] 
```
@inproceedings{roebbenComparative2025,
  author={Roebben, Arnout and van Waterschoot, Toon and Moonen, Marc},
  booktitle={Proc. 2025 IEEE Int. Conf. Acoust., Speech, Signal Process. (ICASSP)}, 
  title={A Comparative Analysis of Generalised Echo and Interference Cancelling and Extended Multichannel Wiener Filtering for Combined Noise Reduction and Acoustic Echo Cancellation}, 
  year={2025},
  keywords={Acoustic echo cancellation (AEC), Noise reduction
(NR), Extended multichannel Wiener filter (MWFext),
Generalised echo and interference canceller (GEIC)}
  }
```

[2]
```
@misc{roebben2025comparativeanalysisgeneralisedecho,
      title={A Comparative Analysis of Generalised Echo and Interference Cancelling and Extended Multichannel Wiener Filtering for Combined Noise Reduction and Acoustic Echo Cancellation}, 
      author={Arnout Roebben and Toon van Waterschoot and Marc Moonen},
      year={2025},
      eprint={2503.03593},
      archivePrefix={arXiv},
      primaryClass={eess.AS},
      url={https://arxiv.org/abs/2503.03593}, 
}
```

[3]
```
@misc{roebbenGithubRepositoryComparative2024,
  title = {Github Repository: {{A Comparative}} Analysis of Generalised Echo and Interference Cancelling and Extended Multichannel {{Wiener}} Filtering for Combined Noise Reduction and Acoustic Echo Cancellation},
  author = {Roebben, A.},
  year = {2024},
  journal = {GitHub},
  urldate = {2023-07-15},
  howpublished = {https://https://github.com/Arnout-Roebben/AECNR\_GEIC\_vs\_MWFext},
  langid = {english}
}
```

[4]
```
@inproceedings{herbordtJointOptimizationLCMV2004,
  title = {Joint Optimization of {{LCMV}} Beamforming and Acoustic Echo Cancellation},
  booktitle = {Proc. 12th  European Signal Process. Conf.  (EUSIPCO)},
  author = {Herbordt, W. and Kellermann, W. and Nakamura, S.},
  year = {2004},
  month = sep,
  pages = {2003--2006},
  address = {Vienna, Austria},
  }
```

[5]
```
@article{ruizDistributedCombinedAcoustic2022,
  title = {Distributed Combined Acoustic Echo Cancellation and Noise Reduction in Wireless Acoustic Sensor and Actuator Networks},
  author = {Ruiz, S. and {van Waterschoot}, T. and Moonen, M.},
  year = {2022},
  journal = {IEEE/ACM Trans. Audio, Speech, Language Process.},
  volume = {30},
  pages = {534--547},
  issn = {2329-9304},
}
```

[6]
```
@inproceedings{markovich-golanPerformanceAnalysisCovariance2015,
  title = {Performance Analysis of the Covariance Subtraction Method for Relative Transfer Function Estimation and Comparison to the Covariance Whitening Method},
  booktitle = {Proc. 2015 IEEE Int. Conf. Acoust., Speech, Signal Process. (ICASSP)},
  author = {{Markovich-Golan}, Shmulik and Gannot, Sharon},
  year = {2015},
  month = apr,
  address = {Brisbane, Australia},
  pages = {544--548},
  issn = {2379-190X},
  doi = {10.1109/ICASSP.2015.7178028},
}
```