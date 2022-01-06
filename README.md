# tsNetworks


Two R script files, containing the procedures to fit non-stationary autoregressive link community models to high-resolution time series of network snapshots. The fitted models provide detailed daily contact probability patterns between actors in the network (also nonstationary ones), and estimates of their correlation (the tendency of an interaction to persist once it is formed) over time. Potential uses include, for example, detailed temporal risk analysis of mobile phone or energy networks for overload, or the spread of an epidemic over a social network.

*EM_BALARM.R*: functions to fit the BALARM (block-autoregressive model) by the EM algorithm.

*Misc.R*: miscellaneous functions for diagnostics and calculations of derived quantities from the fitted models. 

#### Reference:
- Süveges, M. and Olhede, S. Networks with Correlated Edge Processes. Submitted to the Journal of the Royal Statistical Society A, Special Issue on Networks and Society, Dec. 2021.
