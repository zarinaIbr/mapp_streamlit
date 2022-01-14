import pandas as pd
import streamlit as st
import pickle
import base64
import pygsheets
from datetime import datetime
import pytz
from collections import defaultdict
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
import io
import plotly.graph_objs as go

st.markdown(
    """<style>
    [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {width: 370px;}
    [data-testid="stSidebar"][aria-expanded="false"] > div:first-child {width: 1000px; margin-left: -500px;}
    </style>""", unsafe_allow_html=True)

st.markdown("<h2 style='text-align: center; color: black;'>Генерация правил для исправления ААО</h2>",
            unsafe_allow_html=True)

@st.cache(allow_output_mutation=True)
def load_data():
    with open('data_mapp_web_m50.pickle', 'rb') as f:
        data = pickle.load(f)
        return data, data.loc[:, ['query', 'freq']]

def svg_html(reac):
    svg = reac.depict()
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<img src="data:image/svg+xml;base64,%s"/>' % b64
    return html

def update_reaction(old_r, f):
    new_r = old_r.copy()
    for mol in new_r.products:
        mol.remap(f)
    return new_r

def gen_fix_r(reaction, num_bad, num_good):
    fix = {t[0]: t[1] for t in zip(num_bad, num_good)}
    reaction_new = update_reaction(reaction, fix)
    return fix, reaction_new

def load_remapping_rules(reactions):
    rules = []
    for bad, good in reactions:
        if str(bad) != str(good):
            raise ValueError('bad and good reaction should be equal')

        cgr_good, cgr_bad = ~good, ~bad
        gc = cgr_good.augmented_substructure(cgr_good.center_atoms, deep=1)
        bc = cgr_bad.augmented_substructure(cgr_bad.center_atoms, deep=1)
        atoms = set(bc.atoms_numbers + gc.atoms_numbers)

        re_g, re_b, pr_g, pr_b = set(), set(), set(), set()
        for pr in good.reactants:
            re_g.update(pr)
        for pr in good.products:
            pr_g.update(pr)
        for pr in bad.products:
            pr_b.update(pr)
        for pr in bad.reactants:
            re_b.update(pr)
        atoms.update((re_b.difference(pr_b)).intersection(pr_g))
        strange_atoms = pr_b.difference(pr_g)
        atoms.update(strange_atoms)
        bad_query = (cgr_bad).substructure(atoms.intersection(cgr_bad), as_query=True)
        good_query = (cgr_good).substructure(atoms.intersection(cgr_good), as_query=True)

        fix = {}
        for mb, mg in zip(bad.products, good.products):
            fix.update({k: v for k, v in zip(mb, mg) if k != v and k in atoms})  # get fix map
        valid = set(fix).difference(strange_atoms)
        rules.append((bad_query, good_query, fix, valid))
    return rules

def download_to_excel(d_info, type_aam):
    path = 'main_secret.json'
    gc = pygsheets.authorize(service_file=path)
    if type_aam:
        sh = gc.open('good_mapping')
    else:
        sh = gc.open('bad_mapping')
    wks = sh.worksheet_by_title('Лист1')
    df_excel = wks.get_as_df()
    df_new = pd.DataFrame(d_info)
    values = df_new[['rc', 'freq']].values.tolist()
    wks.append_table(values, start=f'A{len(df_excel) + 1}', end=None, dimension='ROWS', overwrite=False)

def download_to_drive(data, time, type_aam):
    with open(f"data_{time}.pickle", 'wb') as f:
        pickle.dump(data, f)
    gauth = GoogleAuth()
    drive = GoogleDrive(gauth)
    upload_file_list = [f"data_{time}.pickle"]
    if type_aam:
        id_folder = '1Yhe2v94GDlqJzNM6K4Gh5Qd2_6tc8Tnd'
    else:
        id_folder = '1HcdJo5VxcpTnVdVT8leG6_bjo3j26CnX'
    for upload_file in upload_file_list:
        gfile = drive.CreateFile({'parents': [{'id': id_folder}]})
        gfile.SetContentFile(upload_file)
        gfile.Upload()

data_all, data_df = load_data()
st.sidebar.write('Датасет', data_df)
selected_indices = st.sidebar.multiselect('Выберите р.ц', data_df.index)
selected_rows = data_df.loc[selected_indices]
st.sidebar.write('Выбранный р.ц:', selected_rows)

st.sidebar.markdown("<h2 style='text-align: center; color: black;'>Инструкция</h2>", unsafe_allow_html=True)
st.sidebar.markdown("""
* 1.Выберите один р.ц по его индексу в таблицe и проверьте его правильность по реакции
* 2.Если маппинг правильный, нажмите Да, и переходите к следующему р.ц
* 3.Если маппинг ошибочный, нажмите Нет
    * 3.1.Если словарь НЕ пустой, нажмите Delete
    * 3.2.Добавляйте пары атомов в пустой словарь, нажимая на Add
* 4.После сбора словаря, проверьте себя по сгенерированной реакции
* 5.Если исправленный маппинг реакции вы считаете верным, нажмите Сгенерировать правило
* 6.После окончания сбора правил, нажмите Save all
""")

if len(selected_indices) == 1:
    reaction = data_all['reaction'].values[int(selected_indices[0])][0]
    st.markdown("<p style='text-align: center; color: black;'>Исходная реакция</p>", unsafe_allow_html=True)
    st.write(svg_html(reaction), unsafe_allow_html=True)

    st.markdown("<p style='text-align: left; color: black;'>ААО верное?</p>", unsafe_allow_html=True)
    checkbox_yes = st.checkbox("Да, корректное")
    checkbox_no = st.checkbox("Нет, ошибочное")
    if checkbox_no and checkbox_yes:
        st.error('Выберите только один вариант')
    elif checkbox_yes and not checkbox_no:
        if 'good_info' not in st.session_state:
            st.session_state.good_info = defaultdict(list)
        tz_NY = pytz.timezone('Europe/Moscow')
        datetime_NY = datetime.now(tz_NY)
        st.session_state.good_info['rc'].append(selected_rows.values.tolist()[0][0])
        st.session_state.good_info['freq'].append(selected_rows.values.tolist()[0][1])
        st.session_state.good_info['nums_reac'].append(data_all['nums_reactions'].values[int(selected_indices[0])])
        time = datetime_NY.strftime("%H_%M")
        b = io.BytesIO()
        pickle.dump(st.session_state.good_info, b)
        b64 = base64.b64encode(b.getvalue()).decode()
        if st.download_button(label="Download data", data=b64, file_name=f"good_data_{time}.pickle"):
            download_to_excel(st.session_state.good_info, type_aam=True)
#         download_to_drive(st.session_state.good_info, datetime_NY.strftime("%H_%M"), type_aam=True)
        del st.session_state.good_info
    elif checkbox_no and not checkbox_yes:
        if 'before' not in st.session_state and 'after' not in st.session_state:
            st.session_state.before, st.session_state.after = [], []
        col1, col2 = st.columns(2)
        with col1:
            sym_b = st.number_input("Ошибочный номер атома", min_value=1)
            sym_a = st.number_input("Корректный номер атома", min_value=1)
        with col2:
            st.write('Сохранить пару')
            add_button = st.button('Add')
            st.write('Очистить словарь')
            del_button = st.button('Delete')
        if add_button:
            st.session_state.before.append(sym_b)
            st.session_state.after.append(sym_a)
        elif del_button:
            while len(st.session_state.after) > 0:
                st.session_state.after.pop()
                st.session_state.before.pop()
        st.code({t[0]: t[1] for t in zip(st.session_state.before, st.session_state.after)})
        col_r, col_rule = st.columns(2)
        with col_r:
            show_reac = st.button('Показать реакцию после исправления ААО')
        with col_rule:
            load_rule = st.button('Сгенерировать правило')
        if show_reac:
            try:
                fix_web, reaction_new = gen_fix_r(reaction, st.session_state.before, st.session_state.after)
                st.write(svg_html(reaction_new), unsafe_allow_html=True)
            except:
                st.error('Ошибка при исправлении ААО')
        if load_rule:
            fix_web, reaction_new = gen_fix_r(reaction, st.session_state.before, st.session_state.after)
            if 'bad_info' not in st.session_state:
                st.session_state.bad_info = defaultdict(list)
            try:
                rule = load_remapping_rules([(reaction, reaction_new)])
                st.session_state.bad_info['rc'].append(selected_rows.values.tolist()[0][0])
                st.session_state.bad_info['freq'].append(selected_rows.values.tolist()[0][1])
                st.session_state.bad_info['nums_reac'].append(data_all['nums_reactions'].values[int(selected_indices[0])])
                st.session_state.bad_info['bad_r'].append(reaction)
                st.session_state.bad_info['good_r'].append(reaction_new)
                st.session_state.bad_info['rule'].append(rule)

                col_bq, col_gq, col_f = st.columns([5, 5, 1])
                with col_bq:                                      
                    st.text('Р.ц непр.ААО')
                    st.write(svg_html(rule[0][0]), unsafe_allow_html=True)
                with col_gq:
                    st.text('Р.ц прав.ААО')
                    st.write(svg_html(rule[0][1]), unsafe_allow_html=True)
                with col_f:
                    st.text('Словарь пар атомов')
                    st.write(f'{rule[0][2]}')
            except:
                st.error('Ошибка при создании правила.')

#         save_all = st.button('Save all')
#         if save_all:
        tz_NY = pytz.timezone('Europe/Moscow')
        datetime_NY = datetime.now(tz_NY)
        time = datetime_NY.strftime("%H_%M")
        b = io.BytesIO()
        pickle.dump(st.session_state.bad_info, b)
        b64 = base64.b64encode(b.getvalue()).decode()
        if st.download_button(label="Download data", data=b64, file_name=f"bad_data_{time}.pickle"):
            download_to_excel(st.session_state.bad_info, type_aam=False)
#             download_to_drive(st.session_state.bad_info, datetime_NY.strftime("%H_%M"), type_aam=False)
        del st.session_state.bad_info

stat = st.button('Посмотреть статистику')
if stat:
    with st.expander("Свернуть/Развернуть", True):
        path = 'main_secret.json'
        gc = pygsheets.authorize(service_file=path)
        sh_good = gc.open('good_mapping')
        sh_bad = gc.open('bad_mapping')
        wks_good = sh_good.worksheet_by_title('Лист1')
        wks_bad = sh_bad.worksheet_by_title('Лист1')
        df_bad, df_good = wks_bad.get_as_df(), wks_good.get_as_df()

        col1, col2, col3 = st.columns(3)
        col1.metric("Всего реакций", "1358000")
        col2.metric("Кол-во правильных реакций", f"{df_good['freq'].sum()}")
        col3.metric("Создано руллов", f"{len(df_bad['freq'])}")

        labels = ['Корректный маппинг', 'Потенциально корректный маппинг', 'Неизвестно']
        values = [df_good['freq'].sum(), df_bad['freq'].sum(), 1358000 - df_good['freq'].sum() - df_bad['freq'].sum()]
        colors = ['green', 'yellow', 'gray']
        fig = go.Figure(data=[go.Pie(labels=labels, values=values)])
        fig.update_traces(marker=dict(colors=colors))
        st.plotly_chart(fig, use_container_width=True)
