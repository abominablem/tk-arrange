# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 19:37:23 2021

@author: marcu
"""

from math import gcd
from functools import reduce
import tkinter as tk

def lcm(*args):
    lcm = 1
    for x in args:
        lcm *= x//gcd(lcm, x)
    return lcm

class WidgetLayout:
    """
    Base class containing layout logic for the WidgetSet classes. Should
    not be instantiated directly.
    """
    def __init__(self, layout):
        self.original_layout = layout[:]
        self.span, self.height, self.indices = None, None, None

        self.layout = self.clean_layout(layout)
        self.indices = self._get_indices(self.layout)
        self.span = self.get_span(self.layout)
        self.height = self.get_height(self.layout)
        self.layout = self.normalise_layout(self.layout)

    def clean_layout(self, layout):
        """
        Clean the specified layout. If all values are non-lists, assume all
        buttons should be placed on a single line. Otherwise, take a non-list
        to mean a widget covering an entire line.
        """
        has_list = False
        for k in layout:
            if isinstance(k, list):
                has_list = True
                break

        if has_list:
            for i, k in enumerate(layout):
                if not isinstance(k, list):
                    layout[i] = [k]
        else:
            layout = [layout]

        return layout

    def normalise_layout(self, layout):
        """
        Normalise the layout so that all lists have the same length. Where
        an impossible layout has been specified, adjust the scaling to
        approximate the input as closely as possible. The heirarchy here is
        to prioritise the scaling in earlier rows and then allocate/trim any
        extra/missing space from re-distributing later button spacing.
        """
        span_dict = {}
        for i, layer in enumerate(layout):
            layer_out = []

            # how much of the width of the layer each item in the list
            # corresponds to, as a fraction
            layer_above = layout[i-1] if i > 0 else []

            for j, x in enumerate(layer):
                if x in layer_out and x > 0:
                    continue

                exp_start = len(layer_out)

                # Spacers have no restrictions on alignment
                if x < 0:
                    act_start = exp_start
                else:
                    try:
                        act_start = layer_above.index(x)
                    except ValueError:
                        act_start = exp_start

                if act_start < exp_start:
                    # if the space being overlaid is all spacers, continue
                    if all(key < 0  for key in layer[act_start:exp_start]):
                        pass
                    else:
                    # TODO - work out how much space must be taken away from
                    # previous buttons. allocate that removed space so that
                    # proportions are preserved
                        raise ValueError(
                            "Unable to draw button set. Objects would overlap")

                # if the actual start is beyond the expected start, then the
                # gap in between must be padded out. Where possible, stretch
                # the button to the left over to fill the gap. Otherwise,
                # fill it with a spacer.
                elif act_start > exp_start:
                    pad_len = act_start - exp_start
                    if len(layer_out) == 0:
                        pad_val = -1
                    elif layer_out[-1] not in layer_above:
                        pad_val = layer_out[-1]
                    else:
                        pad_val = -1
                    layer_out += [pad_val]*pad_len

                    scale = lcm(pad_len, len(layer))
                    for k in range(i):
                        layout[k] = self._scale_list(layout[k], scale)
                    for key in span_dict:
                        span_dict[key] *= scale
                    layer_out = self._scale_list(layer_out, scale)
                    layer_above = self._scale_list(layer_above, scale)
                    self.span *= scale

                available_span = self.span - len(layer_out)
                if x < 0:
                    count = self._count_sequence(layer, x, start = j)
                else:
                    count = layer.count(x)

                exp_span = int(count/len(layer[j:]) * available_span)
                act_span = (span_dict[x] if x in span_dict.keys() and x >= 0
                            else exp_span)
                span_diff = exp_span - act_span

                # allocate the difference if the expected and actual spans
                # are different
                scale = 1
                if span_diff != 0:
                    unallocated_items = [k for k in layer
                                         if not (k in layer_out or k == x)]
                    # scale up the whole layout to ensure the difference can
                    # be properly allocated according to the defined
                    # proportions
                    scale = len(unallocated_items)
                    for j in range(i):
                        layout[j] = self._scale_list(layout[j], scale)
                    for key in span_dict:
                        span_dict[key] *= scale
                    layer_out = self._scale_list(layer_out, scale)
                    layer_above = self._scale_list(layer_above, scale)
                    self.span *= scale

                span_dict[x] = act_span*scale
                layer_out += [x]*(act_span*scale)

            layout[i] = layer_out

        # pad the layout to ensure each layer is the same length
        len_dict = {i: len(layer) for i, layer in enumerate(layout)}
        max_layer_len = max(list(len_dict.values()))
        for i, layer in enumerate(layout):
            layout[i] = layer + [-1]*(max_layer_len - len_dict[i])

        # reduce the layout to the smallest version preserving the allocated
        # spacing
        layout, span_gcd = self._reduce_layout(layout)
        # self.span_dict = {k: int(span_dict[k]/span_gcd) for k in span_dict}
        return layout

    def _scale_list(self, lst, scale):
        lst_out = []
        for x in lst:
            lst_out += [x]*scale
        return lst_out

    def _reduce_layout(self, layout):
        """
        Reduce a layout to the smallest possible version which preserves all
        scaling.
        """
        counts = []
        for i in self._get_indices(layout):
            for layer in layout:
                counts.append(layer.count(i))

        count_gcd = reduce(gcd, counts)
        return [layer[::count_gcd] for layer in layout], count_gcd

    def _count_sequence(self, lst, value, start = 0):
        """
        Count the number of sequential occurrences of a value in a list after
        a defined point.
        """
        lst = lst[start:]
        try:
            lst = lst[lst.index(value):]
            i = 0
            while i < len(lst) and lst[i] == value:
                i += 1
            return i
        except ValueError:
            return 0

    def _get_indices(self, layout):
        if self.indices is None:
            indices = [i for layer in layout for i in layer if i >= 0]
            self.indices = list(set(indices))
        return self.indices

    def _check_layout(self, layout):
        """ for debug purposes only """
        for i, layer in enumerate(layout):
            vdict = {x: layer.count(x) for x in set(layer)}
            print("Layer %s: %s, " % (i, len(layer)), vdict, layer)

        for num in self.indices:
            print({"num": num,
                   "row": self._get_grid_row(num),
                   "column": self._get_grid_column(num),
                   "rowspan": self._get_row_span(num),
                   "columnspan": self._get_column_span(num)
                   })

    def get_span(self, layout = None):
        if self.span is None:
            if layout is None: layout = self.layout
            lengths = [len(k) for k in layout]
            self.span = lcm(*lengths)
        return self.span

    def get_height(self, layout = None):
        if self.height is None:
            if layout is None: layout = self.layout
            self.height = len(layout)
        return self.height



class WidgetSetComponent:
    """
    Container for attributes of a widget. Grid references must be set
    externally.
    """
    def __init__(self, id, wdict):
        self.name = __class__.__name__
        self.id = id
        self.widget = wdict.get("widget", None)
        self.grid_kwargs = wdict.get("grid_kwargs", {})
        self.widget_kwargs = wdict.get("widget_kwargs", {})
        self.placed = False
        self.master = None if self.widget is None else self.widget.master

        self.stretch_height = wdict.get("stretch_height", False)
        self.stretch_width = wdict.get("stretch_width", False)

        self.stretch_height_weight = wdict.get(
            "stretch_height_weight",
            1 if self.stretch_height else 0
            )
        self.stretch_width_weight = wdict.get(
            "stretch_width_weight",
            1 if self.stretch_width else 0
            )

        self.row = None
        self.column = None
        self.rowspan = None
        self.columnspan = None

        self.is_spacer = (self.id < 0)

    def __iter__(self):
        """ For simplicity, create an iterator of length 1 """
        self._i = -1
        return self

    def __next__(self):
        if self._i < 0:
            self._i += 1
            return self
        else:
            raise StopIteration

    def set_grid_refs(self, **kwargs):
        self.__dict__.update(kwargs)

    def grid_refs(self):
        return {"row": self.row,
                "column": self.column,
                "rowspan": self.rowspan,
                "columnspan": self.columnspan}

    def grid(self, **kwargs):
        self.widget.grid(**kwargs)
        self.placed = True

    def check_collision(self, rc):
        assert len(rc) == 2
        # If widget has not been placed in window, it cannot collide
        if not self.placed: return False
        r, c = rc
        row_collides, col_collides = True, True
        if r is not None:
            row_collides = (self.row <= r < self.row + self.rowspan)
        if c is not None:
            col_collides = (self.column <= c < self.column + self.columnspan)
        return row_collides and col_collides


class WidgetSet(WidgetLayout):
    def __init__(self, frame, widgets, layout):
        """
        Parameters
        ----------
        frame : tk.Frame
            frame to arrange widgets inside.
        widgets : dict
            dict of dictionaries with keys the corresponding numbers used
            when specifying the layout
            dict can have keys:
                widget: tk.Widget object
                grid_kwargs: kwargs to pass to tk.Widget.grid method
                widget_kwargs: if widget is created at runtime, kwargs to pass
                               to the widget on creation
                stretch_height: Boolean, allow widget to stretch vertically.
                                Default is False.
                stretch_height_weight: Integer, weighting for vertical space
                                       allocation
                stretch_width: Boolean, allow widget to stretch horizontally.
                               Default is False.
                stretch_width_weight: Integer, weighting for horizontal space
                                       allocation
        layout : list
            Up to 2-dimensional list defining the layout of the buttons. Use
            sequential numbers aligning to the list of button names.
            [[1, 2 , 2],[1, 3, 4],[5]] will result in a set of buttons that
            looks like the following:
            [1][  2 ]
            ["][3][4]
            [   5   ]
            The number of entries in each row will weight the relative width of
            each button.
            Use negative numbers to insert a blank space in the button set.
            Where buttons span multiple rows, the button width is determined by
            the first row it appears.
        frm_kwargs : dict
            kwargs to pass to the frame object.
        """
        super().__init__(layout)
        self.frame = frame
        self.master = frame
        self.wdict = widgets
        self.widgets = {key: WidgetSetComponent(key, widgets[key])
                        for key in widgets if key in self.indices}

        self.create_widgets()

    def create_widgets(self):
        if self.widgets == {}: return
        for wdgt_lst in self.widgets.values():
            for wdgt in wdgt_lst:
                if wdgt.is_spacer: continue
                wdgt.set_grid_refs(**self._get_grid_refs(wdgt.id))
                wdgt.grid(**wdgt.grid_refs(), **wdgt.grid_kwargs)

        self._create_spacers()
        self.rc_configure()

    def _get_grid_refs(self, key):
        return {"row": self._get_grid_row(key),
                "column": self._get_grid_column(key),
                "rowspan": self._get_row_span(key),
                "columnspan": self._get_column_span(key)}

    def _get_row_span(self, key):
        i = 0
        for n, layer in enumerate(self.layout):
            if key in layer:
                i += 1
        return None if i == 0 else i

    def _get_column_span(self, key):
        for layer in self.layout:
            count = layer.count(key)
            if count != 0: return count

    def _get_grid_row(self, key):
        for n, layer in enumerate(self.layout):
            if key in layer:
                return n

    def _get_grid_column(self, key):
        """
        Return index of first column widget appears
        """
        for layer in self.layout:
            try:
                return layer.index(key)
            except:
                continue

    def _create_spacers(self):
        for i, layer in enumerate(self.layout):
            for j, x in enumerate(layer):
                if x < 0:
                    widget = self._get_spacer(x, self.wdict.get(x, {}))
                    widget.set_grid_refs(row = i, column = j,
                                         rowspan = 1, columnspan = 1)
                    widget.grid(**widget.grid_refs(), **widget.grid_kwargs)
                    if not x in self.widgets:
                        self.widgets[x] = [widget]
                    else:
                        self.widgets[x] += widget

    def _get_spacer(self, id, wdict):
        if not 'widget' in wdict:
            widget_kwargs = wdict.get('widget_kwargs', {})
            wdict['widget'] = tk.Label(self.frame, **widget_kwargs)
        return WidgetSetComponent(id, wdict)

    def rc_configure(self):
        """
        Configure the row/column allocation of scaling changes. Strictly
        enforces fixed width assumptions. If unspecified, assume that the
        widget should be fixed width.
        """
        # dictionary of whether to allow each row or column to stretch
        rc_cfg = {
            "column": {c: True for c in range(self.get_span())},
            "row": {r: True for r in range(self.get_height())}
            }

        # Check each cell and fix width/height of those intersecting a fixed
        # width/height widget
        for r, layer in enumerate(self.layout):
            for c, key in enumerate(layer):
                if key < 0: continue

                if not self.widgets[key].stretch_height:
                    rc_cfg["row"][r] = False

                if not self.widgets[key].stretch_width:
                    rc_cfg["column"][c] = False

        for c, bln in rc_cfg["column"].items():
            if not bln:
                self.frame.columnconfigure(c, weight = 0)
            else:
                weights = [0]
                for widget_lst in self.widgets.values():
                    for widget in widget_lst:
                        if widget.check_collision((None, c)):
                            weights.append(widget.stretch_width_weight)
                self.frame.columnconfigure(c, weight = max(weights))

        for r, bln in rc_cfg["row"].items():
            if not bln:
                self.frame.rowconfigure(r, weight = 0)
            else:
                weights = [0]
                for widget_lst in self.widgets.values():
                    for widget in widget_lst:
                        if widget.check_collision((r, None)):
                            weights.append(widget.stretch_height_weight)
                self.frame.rowconfigure(r, weight = max(weights))

        self.rc_config = rc_cfg

    def grid(self, **kwargs):
        self.frame.grid(**kwargs)

    def rowconfigure(self, *args, **kwargs):
        self.frame.rowconfigure(*args, **kwargs)

    def columnconfigure(self, *args, **kwargs):
        self.frame.columnconfigure(*args, **kwargs)


class ButtonSet(WidgetSet):
    """
    Generate a bound button set from a dictionary of buttons and bindings and
    corresponding layout matrix.

    Creates a frame which contains the button set.
    """
    def __init__(self, master, buttons, layout, frm_kwargs, set_width = None):
        self.master = master
        self.buttons = buttons
        self.widgets = {}
        self.set_width = set_width
        self.frame = tk.Frame(master, **frm_kwargs)
        super().__init__(frame = self.frame, widgets = {}, layout = layout)

        for key in self._get_indices(self.layout):
            if key < 0:
                if not self.set_width is None:
                    self.buttons[key] = {
                        "widget_kwargs": {
                            "width": int(self.set_width/self.span)
                            }
                        }
            else:
                self.add_button(key)
        self.set_widgets(self.buttons)
        self.create_widgets()

    def set_widgets(self, wdict):
        self.widgets = {key: WidgetSetComponent(key, wdict[key])
                        for key in wdict}

    def add_button(self, key):
        if key < 0: return

        btn_args = self.buttons[key]
        btn_args.setdefault("widget_kwargs", {})
        if not self.set_width is None:
            btn_width = int(self.set_width*self._get_column_span(key)/self.span)
            btn_args["widget_kwargs"]["width"] = btn_width

        # Create tk button
        btn = tk.Button(
            master = self.frame,
            text = btn_args["label"],
            **btn_args["widget_kwargs"]
            )

        bindings = btn_args.get("bindings")
        if bindings is not None:
            events = bindings["event"]
            functions = bindings["function"]
            if not isinstance(events, list): events = [events]
            if not isinstance(functions, list): functions = [functions]
            for i in range(len(events)):
                event = events[i]
                function = functions[i]
                btn.bind(event, function)

        self.buttons[key]["widget"] = btn
        try:
            self.__dict__[self.buttons[key]["name"]] = btn
        except KeyError:
            pass