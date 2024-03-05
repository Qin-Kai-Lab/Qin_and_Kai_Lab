import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import os

df = pd.read_csv("RNA_CellsPerClust.csv")

targDir = './Paper_figs/Fig1/Cell_Numbers/'

fig = go.Figure(data=[go.Pie(labels=df["Cluster"], values=df["Freq"])])
fig.update_layout(title='Pie Chart')

fig.update_layout(
    title={
        'text': 'Number of Cells Per Cluster',
        'x': 0.5,  # Center align the title
        'y': 0.95  # Position the title at the top
    },
    width=1000,  # Set the width of the figure to 1000 pixels
    height=1000,
    legend=dict(
        font=dict(
            size=20  # Set the font size of the legend
        )
    )
)

# Save the figure with specified dimensions and DPI
pio.write_image(fig, f'{targDir}RNA_CellsPerClust_Pie.png', width=1000, height=1000, scale=3)
