from wordcloud import WordCloud, get_single_color_func
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


class SimpleGroupedColorFunc(object):
    """Create a color function object which assigns EXACT colors
       to certain words based on the color to words mapping

       Parameters
       ----------
       color_to_words : dict(str -> list(str))
         A dictionary that maps a color to the list of words.

       default_color : str
         Color that will be assigned to a word that's not a member
         of any value from color_to_words.
    """

    def __init__(self, color_to_words, default_color="grey"):
        self.word_to_color = {word: color
                              for (color, words) in color_to_words.items()
                              for word in words}

        self.default_color = default_color

    def __call__(self, word, **kwargs):
        return self.word_to_color.get(word, self.default_color)


class GroupedColorFunc(object):
    """Create a color function object which assigns DIFFERENT SHADES of
       specified colors to certain words based on the color to words mapping.

       Uses wordcloud.get_single_color_func

       Parameters
       ----------
       color_to_words : dict(str -> list(str))
         A dictionary that maps a color to the list of words.

       default_color : str
         Color that will be assigned to a word that's not a member
         of any value from color_to_words.
    """

    def __init__(self, color_to_words, default_color = "grey"):
        self.color_func_to_words = [
            (get_single_color_func(color), set(words))
            for (color, words) in color_to_words.items()]

        self.default_color_func = get_single_color_func(default_color)

    def get_color_func(self, word):
        """Returns a single_color_func associated with the word"""
        try:
            color_func = next(
                color_func for (color_func, words) in self.color_func_to_words
                if word in words)
        except StopIteration:
            color_func = self.default_color_func

        return color_func

    def __call__(self, word, **kwargs):
        return self.get_color_func(word)(word, **kwargs)


def draw_word_cloud(df, score='OddRatio', N=20, N_under_represented=0, under_represented=False, scale=2.0):
    """
    Draw first and last N descriptions of df
    Scale less abboundant terms by a factor scale
    """
    if under_represented:
        # draw also underepresented words
        head, tail = df[:N].copy(), df[-N_under_represented:].copy()
        
        # invert score of tail dataframe
        if score=='OddRatio':
            tail[score] = 1 / tail[score] * scale
        
        df_draw = pd.concat([head,tail], ignore_index=True)
    else:
        df_draw = df[:N]
    
    wc = WordCloud(background_color='white', width=1600, height=800, max_words=100, collocations=False)
    
    if score=='p-value':
        wc.generate_from_frequencies(
            {line['label']: np.log(1/line[score]) for _, line in df_draw.iterrows()}
        )
    else:
        wc.generate_from_frequencies(
            {line['label']: line[score] for _, line in df_draw.iterrows()}
        )

    if under_represented:
        color_to_words = {
            '#ff7f0e' : head['label'].to_list(),
            '#1f77b4' : tail['label'].to_list()
        }
        # Create a color function with single tone
        grouped_color_func = SimpleGroupedColorFunc(color_to_words)

        # Apply our color function
        wc.recolor(color_func=grouped_color_func)

        return wc
    
    else:
        return wc
    